use crate::io::FqReader;
use crate::qual_map::QualMap;
use crate::{dna, io, stats};

use bio::io::fastq::Record;
use ndarray::{Array2, Axis};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn fq_check(
    fq_path: &str,
    out_path: &str,
    qual_threshold: u8,
    ascii_bases: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    if ascii_bases != 33u8 && ascii_bases != 64u8 {
        println!("Warning: Input ascii base {} is not 33 or 64.", ascii_bases);
    }
    let mut out_path: String = out_path.to_string();
    out_path.push_str(".tsv");

    let (max_len, qual_idx_map) =
        parse_fq_and_write_length_and_qual_stats(fq_path, &out_path, qual_threshold, ascii_bases)?;

    process_fastq_and_check(
        fq_path,
        &out_path,
        qual_threshold,
        ascii_bases,
        max_len,
        &qual_idx_map,
    )?;
    Ok(())
}

fn process_fastq_and_check(
    fq_path: &str,
    out_path: &str,
    qual_threshold: u8,
    ascii_bases: u8,
    max_length: usize,
    qual_map: &HashMap<u8, usize>,
) -> Result<(), std::io::Error> {
    let mut dna_count_arr = Array2::<f64>::zeros((max_length, 5));
    let mut qual_sum_arr = Array2::<f64>::zeros((max_length, 2));
    let mut qual_count_arr =
        QualMap::init_qual_count_arr(max_length, qual_map.len(), qual_threshold);

    let fq_iter = FqReader::new(fq_path)?;
    for record in fq_iter.records() {
        let cur_read = record.unwrap();

        update_dna_count(&mut dna_count_arr, &cur_read);
        update_qual_sum(&mut qual_sum_arr, &cur_read);
        qual_count_arr.update_qual_count_mat(&cur_read, qual_map, qual_threshold + ascii_bases);
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
) -> Result<(), std::io::Error> {
    let mut writer = io::append_bufwriter(path)?;
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
fn parse_fq_and_write_length_and_qual_stats(
    fq_path: &str,
    out_path: &str,
    qual_threshold: u8,
    ascii_bases: u8,
) -> Result<(usize, HashMap<u8, usize>), Box<dyn std::error::Error>> {
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
        qual_set: HashSet<u8>,
        ascii_bases: u8,
    ) -> (HashMap<u8, usize>, String) {
        let mut colname = format!("POS\t#bases{}\tavgQ\terrQ", dna::ACGTN_TAB);
        let mut qual_idx_map: HashMap<u8, usize> = HashMap::new();
        if qual_threshold == 0 {
            let mut qual_vec: Vec<u8> = qual_set.into_iter().collect();
            qual_vec.sort();
            let mut ordered_qual = String::new();
            for (idx, &qual) in qual_vec.iter().enumerate() {
                qual_idx_map.insert(qual, idx);
                ordered_qual.push_str(&format!("\t%{}", qual - ascii_bases));
            }
            colname.push_str(&ordered_qual);
        } else {
            colname.push_str("\t%low\t%high");
        }
        colname.push('\n');
        (qual_idx_map, colname)
    }
    let mut reads_size = Vec::<u32>::new();
    let mut qual_set: HashSet<u8> = HashSet::new();
    let fq_iter = FqReader::new(fq_path)?;
    for record in fq_iter.records() {
        let cur_read = record.unwrap();
        reads_size.push(cur_read.seq().len() as u32);
        qual_set.extend(cur_read.qual());
    }
    reads_size.par_sort_unstable();
    let mut out_string = get_read_len_stats_str(&reads_size);

    out_string.push_str(&format!(
        "[Quality value] {} distinct quality values.\n",
        qual_set.len()
    ));

    let (qual_idx_map, colname) = get_colname_str(qual_threshold, qual_set, ascii_bases);
    out_string.push_str(&colname);

    let file = File::create(out_path)?;
    let mut buf_writer = BufWriter::new(file);
    write!(&mut buf_writer, "{}", out_string)?;

    Ok((*reads_size.last().unwrap() as usize, qual_idx_map))
}
fn update_dna_count(dna_count_mat: &mut Array2<f64>, read: &Record) {
    read.seq().iter().enumerate().for_each(|(pos, &c)| {
        let char_idx = dna::get_dna_idx_from_u8(c);
        dna_count_mat[(pos, char_idx)] += 1.0;
    });
}
fn update_qual_sum(qual_mat: &mut Array2<f64>, read: &Record) {
    for (i, &q) in read.qual().iter().enumerate() {
        let q_f64 = q as f64;
        qual_mat[(i, 0)] += stats::convert_q_score_to_p_err(q_f64);
        qual_mat[(i, 1)] += q_f64;
    }
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
    let p_avg = qual_sum_mat[(upos, 0)] / total;
    let q_avg = qual_sum_mat[(upos, 1)] / total;
    let q_err = stats::convert_p_err_to_q_score(p_avg);
    format!("\t{:.1}\t{:.1}", q_avg, q_err)
}
/// Get {avgQ}\t{Qerr} String from qual_sum_map
fn get_total_qual_sum_stats(qual_sum_mat: &Array2<f64>, total: f64) -> String {
    let total_qual_sum = qual_sum_mat.sum_axis(Axis(0)); // Sum of Q/P f across rows
    let p_avg = total_qual_sum[0] / total;
    let q_avg = total_qual_sum[1] / total;
    let q_err: f64 = stats::convert_p_err_to_q_score(p_avg);
    format!("\t{:.1}\t{:.1}", q_avg, q_err)
}
