use crate::dna::DNA;
use crate::fx_iterator::FxIterator;
use crate::qual_map::QualMap;
use crate::sequence_read::SequenceRead;
use crate::stats;

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, prelude::*, BufWriter};

pub fn fq_check(fq_path: &str, qual_threshold: u8, ascii_bases: u8) {
    if ascii_bases != 33u8 && ascii_bases != 64u8 {
        println!("Warning: Input ascii base {} is not 33 or 64.", ascii_bases);
    }
    process_fastq_and_check(fq_path, qual_threshold, ascii_bases)
        .expect("Failed to process fastq file.");
}

fn process_fastq_and_check(fq_path: &str, qual_threshold: u8, ascii_bases: u8) -> io::Result<()> {
    let mut reads_size = Vec::<u32>::new();
    let mut qual_value_set: HashSet<u8> = HashSet::new();
    let mut dna_count_map = QualMap::init_dna_count_map();
    let mut qual_sum_map = QualMap::init_qual_sum_map();
    let mut qual_count_map = QualMap::init_qual_count_map(qual_threshold);

    let fq_iter = FxIterator::new(fq_path, "fastq").unwrap();
    for fx_lines in fq_iter {
        let cur_read = create_sequenceread(fx_lines, ascii_bases)?;
        reads_size.push(cur_read.get_seq_length() as u32);
        qual_value_set.extend(cur_read.get_q_score_vec());
        QualMap::update_map_with_read(&mut dna_count_map, &cur_read, qual_threshold);
        QualMap::update_map_with_read(&mut qual_sum_map, &cur_read, qual_threshold);
        QualMap::update_map_with_read(&mut qual_count_map, &cur_read, qual_threshold);
    }

    write_fq_stat_to_file(
        "test.txt",
        &dna_count_map,
        &qual_count_map,
        &qual_sum_map,
        &reads_size,
        qual_value_set.len(),
        qual_threshold,
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
    dna_count_map: &QualMap,
    qual_count_map: &QualMap,
    qual_sum_map: &QualMap,
    reads_size: &[u32],
    uniq_qual_value_count: usize,
    qual_threshold: u8,
) -> io::Result<()> {
    let col_name = get_colname(qual_threshold);
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "[Quality value] {} distinct quality values.",
        uniq_qual_value_count
    )?;
    writeln!(writer, "{}", get_read_len_stats_string(reads_size))?;
    writeln!(writer, "{}", col_name)?;
    // for pos in 0..(dna_count_map.len()) {
    //     writeln!(
    //         writer,
    //         "{}",
    //         get_dna_stats_string(
    //             dna_count_map,
    //             qual_count_map,
    //             qual_sum_map,
    //             pos
    //         )
    //     )?;
    // }
    Ok(())
}
fn get_colname(qual_threshold: u8) -> String {
    fn get_ordered_dna() -> String {
        let mut out = String::new();
        for dna in DNA::all_chars() {
            out.push_str(&format!("\t{}", dna));
        }
        out
    }
    let mut colname = String::new();
    colname.push_str("POS\t#bases");
    colname.push_str(&get_ordered_dna());
    colname.push_str("\tavgQ\terrQ");
    if qual_threshold == 0 {
        colname.push_str("");
    } else {
        colname.push_str("\t%low\t%high");
        colname.push_str(&get_ordered_dna());
    }
    colname
}
fn get_read_len_stats_string(reads_size: &[u32]) -> String {
    let mut sorted_reads_size = reads_size.to_vec();
    sorted_reads_size.sort();
    format!(
        "[Length] mean={} ; min={} ; med={} ; max={} ; N50={}",
        stats::average(&sorted_reads_size),
        sorted_reads_size.first().unwrap(),
        stats::median_sorted_arr(&sorted_reads_size),
        sorted_reads_size.last().unwrap(),
        stats::n50_sorted_arr(&sorted_reads_size).unwrap()
    )
}
fn get_dna_stats_string(
    dna_count_map: &QualMap,
    qual_sum_map: &QualMap,
    qual_count_map: &QualMap,
    pos: usize,
) -> String {
    let total = QualMap::sum_map_count(dna_count_map, pos) as f64;

    let output_dna_stats = QualMap::get_dna_count_stats(dna_count_map, pos, total);
    let qual_sum_stats = QualMap::get_qual_sum_stats(qual_sum_map, pos, total);
    // let low_qual_count_stats = QualMap::get_low_qual_stats(low_qual_count_map,pos, total);
    let mut out = String::new();
    out.push_str(&output_dna_stats);
    out.push_str(&qual_sum_stats);
    out
}
