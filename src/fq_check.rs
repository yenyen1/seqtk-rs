use crate::dna::DNA;
use crate::fx_iterator::FxIterator;
use crate::qual_map::QualMap;
use crate::sequence_read::SequenceRead;
use crate::stats;

use ndarray::{Array2, Axis};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{self, prelude::*, BufWriter};
use std::time::Instant;

pub fn fq_check(fq_path: &str, qual_threshold: u8, ascii_bases: u8) -> io::Result<()> {
    let out_path = "test.txt";
    let start = Instant::now();

    if ascii_bases != 33u8 && ascii_bases != 64u8 {
        println!("Warning: Input ascii base {} is not 33 or 64.", ascii_bases);
    }

    let (max_len, qual_idx_map) =
        parse_fq_and_write_length_and_qual_stats(fq_path, out_path, qual_threshold, ascii_bases)?;
    process_fastq_and_check(
        fq_path,
        out_path,
        qual_threshold,
        ascii_bases,
        max_len,
        qual_idx_map,
    )?;

    println!(
        "Process time: {:?}",
        format_duration(start.elapsed().as_secs())
    );
    Ok(())
}

fn process_fastq_and_check(
    fq_path: &str,
    out_path: &str,
    qual_threshold: u8,
    ascii_bases: u8,
    max_length: usize,
    qual_map: HashMap<u8, usize>,
) -> io::Result<()> {
    let mut dna_count_arr = Array2::<f64>::zeros((max_length, 5));
    let mut qual_sum_arr = Array2::<f64>::zeros((max_length, 2));
    let mut qual_count_arr =
        QualMap::init_qual_count_arr(max_length, qual_map.len(), qual_threshold);

    let fq_iter = FxIterator::new(fq_path, "fastq").unwrap();
    for fx_lines in fq_iter {
        let cur_read = SequenceRead::create_read(fx_lines[1].trim(), fx_lines[3].trim(),ascii_bases);

        update_dna_count(&mut dna_count_arr, &cur_read); // qual_value=0: count dna without threshold
        update_qual_sum(&mut qual_sum_arr, &cur_read);
        qual_count_arr.update_qual_count_mat(&cur_read, &qual_map, qual_threshold);
    }

    write_fq_stat_to_file(
        out_path,
        &dna_count_arr,
        &qual_count_arr,
        &qual_sum_arr,
        max_length,
    )?;

    Ok(())
}

/// [[[ Write stats files ]]]
fn write_fq_stat_to_file(
    path: &str,
    dna_count_mat: &Array2<f64>,
    qual_count_mat: &QualMap,
    qual_sum_map: &Array2<f64>,
    max_length: usize,
) -> io::Result<()> {
    let mut writer = append_bufwriter(path)?;
    let (total, dna_count_str) = get_total_and_dna_proportion(dna_count_mat);
    let qual_sum_str = get_total_qual_sum_stats(qual_sum_map, total);
    let qual_count_str = qual_count_mat.get_total_qual_count_stats(total);
    writeln!(
        writer,
        "{}{}{}",
        dna_count_str, qual_sum_str, qual_count_str
    )?;
    for pos in 0..max_length {
        let (total, dna_count_str) = get_pos_total_and_dna_proportion(dna_count_mat, pos);
        let qual_sum_str = get_qual_sum_stats(qual_sum_map, pos, total);
        let qual_count_str = qual_count_mat.get_qual_count_stats(pos, total);
        writeln!(
            writer,
            "{}{}{}",
            dna_count_str, qual_sum_str, qual_count_str
        )?;
    }
    Ok(())
}

fn format_duration(duration: u64) -> String {
    let seconds = duration;
    let minutes = seconds / 60;
    let hours = minutes / 60;
    let days = hours / 24;

    let formatted_time = format!(
        "{}-{:02}:{:02}:{:02}",
        days,
        hours % 24,
        minutes % 60,
        seconds % 60
    );

    formatted_time
}
fn append_bufwriter(out_path: &str) -> io::Result<BufWriter<File>> {
    let file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Open the file in append mode (no truncation)
        .open(out_path)?; // Open the file at the provided path
    let writer = BufWriter::new(file);
    Ok(writer)
}

fn parse_fq_and_write_length_and_qual_stats(
    fq_path: &str,
    out_path: &str,
    qual_threshold: u8,
    ascii_bases: u8
) -> Result<(usize, HashMap<u8, usize>), io::Error> {
    fn get_read_len_stats_str(sorted_reads_size: &[u32]) -> String {
        format!(
            "[Length] {} reads; mean={:2} ; min={} ; med={} ; max={} ; N50={}\n",
            sorted_reads_size.len(),
            stats::average(sorted_reads_size),
            sorted_reads_size.first().unwrap(),
            stats::median_sorted_arr(sorted_reads_size),
            sorted_reads_size.last().unwrap(),
            stats::n50_sorted_arr(sorted_reads_size).unwrap()
        )
    }
    fn get_colname_str(
        qual_threshold: u8,
        qual_set: HashSet<char>,
        ascii_bases: u8,
    ) -> (HashMap<u8, usize>, String) {
        let ordered_dna: String = DNA::all_chars().iter().map(|d| format!("\t%{}", d)).fold(
            String::new(),
            |mut acc, s| {
                acc.push_str(&s);
                acc
            },
        );

        let mut colname = format!("POS\t#bases{}\tavgQ\terrQ", ordered_dna);
        let mut qual_idx_map: HashMap<u8, usize> = HashMap::new();
        if qual_threshold == 0 {
            let mut qual_vec: Vec<u8> = qual_set.into_iter().map(|c| c as u8 - ascii_bases).collect();
            qual_vec.sort();
            let mut ordered_qual = String::new();
            for (idx, &qual) in qual_vec.iter().enumerate() {
                qual_idx_map.insert(qual, idx);
                ordered_qual.push_str(&format!("\t%{}", qual));
            }
            colname.push_str(&ordered_qual);
        } else {
            colname.push_str("\t%low\t%high");
        }
        colname.push('\n');
        (qual_idx_map, colname)
    }
    let fq_iter = FxIterator::new(fq_path, "fastq").unwrap();
    let mut reads_size = Vec::<u32>::new();
    let mut qual_set: HashSet<char> = HashSet::new();
    for lines in fq_iter {
        reads_size.push(lines[1].trim().len() as u32);
        qual_set.extend(lines[3].trim().chars());
    }
    reads_size.par_sort_unstable();
    let mut out_string = get_read_len_stats_str(&reads_size);
    out_string.push_str(&format!("[Quality value] {} distinct quality values.\n", qual_set.len()));

    let (qual_idx_map, colname) = get_colname_str(qual_threshold, qual_set, ascii_bases);
    out_string.push_str(&colname);


    let file = File::create(out_path)?;
    let mut buf_writer = BufWriter::new(file);
    write!(&mut buf_writer, "{}", out_string)?;

    Ok((*reads_size.last().unwrap() as usize, qual_idx_map))
}
fn update_dna_count(dna_count_mat: &mut Array2<f64>, read: &SequenceRead) {
    for (pos, &idx) in read.get_seq_idx().iter().enumerate() {
        if let Some(value) = dna_count_mat.get_mut((pos, idx)) {
            *value += 1.0;
        }
    }
}
fn update_qual_sum(qual_mat: &mut Array2<f64>, read: &SequenceRead) {
    fn add(
        qual_mat: &mut Array2<f64>,
        qual_vec: &[f64],
        type_idx: usize, // P=0; Q=1
    ) {
        for (i, &v) in qual_vec.iter().enumerate() {
            if let Some(value) = qual_mat.get_mut((i, type_idx)) {
                *value += v;
            }
        }
    }
    let p_err_vec = read.get_p_err_vec(); // P_error
    add(qual_mat, &p_err_vec, 0);
    let q_score_vec = read.get_q_score_vec_f64(); // Q score
    add(qual_mat, &q_score_vec, 1);
}
/// Get Total\t{#Base}\t{%A}\t{%C}\t{%G}\t{%T}\t{%N} String from dna_count_mat
fn get_total_and_dna_proportion(dna_count_mat: &Array2<f64>) -> (f64, String) {
    let total_dna_counts: Vec<f64> = (dna_count_mat)
        .axis_iter(Axis(1)) // Sum across rows
        .into_par_iter()
        .map(|col| col.sum())
        .collect();
    let total: f64 = total_dna_counts.iter().sum();
    let prop_str: String = total_dna_counts
        .iter()
        .map(|&c| format!("\t{:.1}", 100.0 * c / total))
        .fold(String::new(), |mut acc, s| {
            acc.push_str(&s);
            acc
        });

    (total, format!("Total\t{:.0}{}", total, prop_str))
}
/// Get {Pos}\t{#Base}\t{%A}\t{%C}\t{%G}\t{%T}\t{%N} String from dna_count_mat
fn get_pos_total_and_dna_proportion(dna_count_mat: &Array2<f64>, upos: usize) -> (f64, String) {
    let pos_total: f64 = dna_count_mat.slice(ndarray::s![upos, ..]).iter().sum(); // Sum of counts at position
    let prop_str: String = dna_count_mat
        .slice(ndarray::s![upos, ..])
        .iter()
        .map(|&c| format!("\t{:.1}", 100.0 * c / pos_total))
        .fold(String::new(), |mut acc, s| {
            acc.push_str(&s);
            acc
        });
    (
        pos_total,
        format!("{}\t{}{}", upos + 1, pos_total, prop_str),
    )
}
/// Get {avgQ}\t{Qerr} String from qual_sum_map
fn get_qual_sum_stats(qual_sum_mat: &Array2<f64>, upos: usize, total: f64) -> String {
    let q_avg = qual_sum_mat.get((upos, 1)).unwrap_or(&0.0) / total;
    let p_avg = qual_sum_mat.get((upos, 0)).unwrap_or(&0.0) / total;
    let q_err = SequenceRead::convert_p_err_to_q_score(p_avg);
    format!("\t{:.1}\t{:.1}", q_avg, q_err)
}
/// Get {avgQ}\t{Qerr} String from qual_sum_map
fn get_total_qual_sum_stats(qual_sum_mat: &Array2<f64>, total: f64) -> String {
    let total_qual_sum = qual_sum_mat.sum_axis(Axis(0)); // Sum of Q/P f across rows
    let q_avg = total_qual_sum.get(1).unwrap_or(&0.0) / total;
    let p_avg = total_qual_sum.get(0).unwrap_or(&0.0) / total;
    let q_err: f64 = SequenceRead::convert_p_err_to_q_score(p_avg);
    format!("\t{:.1}\t{:.1}", q_avg, q_err)
}
