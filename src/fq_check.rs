// use rayon::iter::ParallelBridge;
use crate::dna::DNA;
use crate::fx_iterator::FxIterator;
use crate::qual_map::QualMap;
use crate::sequence_read::SequenceRead;
use crate::stats;
use rayon::prelude::*;

// use std::fmt::{Debug, Display};
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
        parse_fq_and_write_length_and_qual_stats(fq_path, out_path, qual_threshold)?;
    process_fastq_and_check(
        fq_path,
        out_path,
        qual_threshold,
        ascii_bases,
        max_len,
        qual_idx_map,
    )
    .expect("Failed to process fastq file.");

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
    qual_set: HashMap<u8, usize>,
) -> io::Result<()> {
    let mut dna_count_mat = QualMap::init_dna_count_mat(max_length);
    let mut qual_sum_map: QualMap = QualMap::init_qual_sum_mat(max_length);
    // let mut qual_count_map = QualMap::init_qual_count_map(qual_threshold);

    let fq_iter = FxIterator::new(fq_path, "fastq").unwrap();
    for fx_lines in fq_iter {
        let cur_read = create_sequenceread(fx_lines, ascii_bases)?;

        QualMap::update_map_with_read(&mut dna_count_mat, &cur_read, qual_threshold); // qual_value=0: count dna without threshold
        QualMap::update_map_with_read(&mut qual_sum_map, &cur_read, qual_threshold);
        // QualMap::update_map_with_read(&mut qual_count_map, &cur_read, qual_threshold);
    }

    write_fq_stat_to_file(
        out_path,
        &dna_count_mat,
        // &qual_count_map,
        &qual_sum_map,
        max_length as usize,
    )?;

    Ok(())
}

/// [[[ Create SequenceRead from fastq read ]]]
fn create_sequenceread(fx_lines: Vec<String>, ascii_bases: u8) -> Result<SequenceRead, io::Error> {
    let mut cur_read = SequenceRead::empty_read(ascii_bases);
    cur_read.set_read_name(fx_lines[0].trim());
    cur_read.set_read_seq(fx_lines[1].trim());
    cur_read.set_read_qual_str(fx_lines[3].trim());
    Ok(cur_read)
}

/// [[[ Write stats files ]]]
fn write_fq_stat_to_file(
    path: &str,
    dna_count_mat: &QualMap,
    // qual_count_map: &QualMap,
    qual_sum_map: &QualMap,
    max_length: usize,
) -> io::Result<()> {
    let mut writer = append_bufwriter(path)?;
    let (total, output_dna_stats) = QualMap::get_total_and_dna_proportion(dna_count_mat);
    for pos in 0..max_length {
        writeln!(
            writer,
            "{}",
            QualMap::get_output(dna_count_mat, qual_sum_map, pos) // get_dna_stats_string(dna_count_map, qual_sum_map, qual_count_map, pos)
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

fn new_bufwriter(out_path: &str) -> io::Result<(BufWriter<File>)> {
    let file = File::create(out_path)?;
    let writer = BufWriter::new(file);
    Ok(writer)
}

fn append_bufwriter(out_path: &str) -> io::Result<(BufWriter<File>)> {
    let file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Open the file in append mode (no truncation)
        .write(true) // We want to write to the file
        .open(out_path)?; // Open the file at the provided path

    let writer = BufWriter::new(file);

    Ok(writer)
}

fn parse_fq_and_write_length_and_qual_stats(
    fq_path: &str,
    out_path: &str,
    qual_threshold: u8,
) -> Result<(usize, HashMap<u8, usize>), io::Error> {
    fn get_read_len_stats_str(sorted_reads_size: &[u32]) -> String {
        format!(
            "[Length] {} reads; mean={:2} ; min={} ; med={} ; max={} ; N50={}\n",
            sorted_reads_size.len(),
            stats::average(&sorted_reads_size),
            sorted_reads_size.first().unwrap(),
            stats::median_sorted_arr(&sorted_reads_size),
            sorted_reads_size.last().unwrap(),
            stats::n50_sorted_arr(&sorted_reads_size).unwrap()
        )
    }
    fn get_uniq_qual_num_str(qual_count: usize) -> String {
        format!("[Quality value] {} distinct quality values.\n", qual_count)
    }
    fn get_colname_str(
        qual_threshold: u8,
        qual_set: HashSet<char>,
    ) -> (HashMap<u8, usize>, String) {
        let ordered_dna: String = DNA::all_chars()
            .iter()
            .map(|d| format!("\t%{}", d))
            .collect();
        let mut colname = format!("POS\t#bases{}\tavgQ\terrQ", ordered_dna);
        let mut qual_idx_map: HashMap<u8, usize> = HashMap::new();
        if qual_threshold == 0 {
            let mut qual_vec: Vec<u8> = qual_set.into_iter().map(|c| c as u8).collect();
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
    out_string.push_str(&get_uniq_qual_num_str(qual_set.len()));

    let (qual_idx_map, colname) = get_colname_str(qual_threshold, qual_set);
    out_string.push_str(&colname);

    let mut buf_writer: BufWriter<File> = new_bufwriter(out_path)?;
    write!(&mut buf_writer, "{}", out_string)?;

    Ok((*reads_size.last().unwrap() as usize, qual_idx_map))
}
