use crate::dna::DNA;
use crate::fx_iterator::FxIterator;
use crate::sequence_read::SequenceRead;
use crate::stats;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, prelude::*};

pub fn fq_check(fq_path: &str, qual_threshold: u8, ascii_bases: u8) {
    if ascii_bases != 33u8 && ascii_bases != 64u8 {
        println!("Warning: Input ascii base {} is not 33 or 64.", ascii_bases);
    }
    process_fastq_and_check(fq_path, qual_threshold, ascii_bases)
        .expect("Failed to process fastq file.");
}

fn process_fastq_and_check(fq_path: &str, qual_threshold: u8, ascii_bases: u8) -> io::Result<()> {
    let mut reads_size = Vec::<u32>::new();
    let mut hashmap_dna_count = init_pos_hashmap_for_dna_count_hashmap(); // 0: Total
    let mut hashmap_low_qual_dna_count = init_pos_hashmap_for_dna_count_hashmap(); // 0: Total
    let mut qual_value_hashset: HashSet<char> = HashSet::new();

    let fq_iter = FxIterator::new(fq_path, "fastq").unwrap();
    for fx_lines in fq_iter {
        let cur_read = create_sequenceread(fx_lines, ascii_bases)?;

        add_dna_count_to_hashmap(&mut hashmap_dna_count, &cur_read);
        add_low_qual_dna_count_to_hashmap(
            &mut hashmap_low_qual_dna_count,
            &cur_read,
            &qual_threshold,
        );
        qual_value_hashset.extend(cur_read.get_qual_chars());
        reads_size.push(cur_read.get_seq_length() as u32);
    }

    write_fq_stat_to_file(
        "test.txt",
        &hashmap_dna_count,
        &hashmap_low_qual_dna_count,
        &reads_size,
        qual_value_hashset.len(),
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

/// [[[ Use hashmap to restore DNA count on each position ]]]
/// Initial and update hashmap
fn add_dna_count_to_hashmap(
    hashmap_dna_count: &mut HashMap<usize, HashMap<DNA, usize>>,
    read: &SequenceRead,
) {
    for (i, char) in read.get_seq().iter().enumerate() {
        add_dna_count_for_a_pos(0, hashmap_dna_count, char); // Total : 0
        add_dna_count_for_a_pos(i + 1, hashmap_dna_count, char);
    }
}
fn add_low_qual_dna_count_to_hashmap(
    hashmap_low_qual_dna_count: &mut HashMap<usize, HashMap<DNA, usize>>,
    read: &SequenceRead,
    qual_threshold: &u8,
) {
    if *qual_threshold != 0_u8 {
        for (i, qual) in read.get_q_qual_score_u8().iter().enumerate() {
            if qual < qual_threshold {
                let dna = read.get_seq_char(i);
                add_dna_count_for_a_pos(0, hashmap_low_qual_dna_count, &dna);
                add_dna_count_for_a_pos(i + 1, hashmap_low_qual_dna_count, &dna);
            }
        }
    }
}
fn init_pos_hashmap_for_dna_count_hashmap() -> HashMap<usize, HashMap<DNA, usize>> {
    // Use pos = 0 to restore total count
    HashMap::from([(0, create_a_dna_to_count_hashmap())])
}
fn create_a_dna_to_count_hashmap() -> HashMap<DNA, usize> {
    let mut dna_hashmap = HashMap::new();
    for c in DNA::all_variants() {
        dna_hashmap.insert(c, 0);
    }
    dna_hashmap
}
fn add_dna_count_for_a_pos(
    pos: usize,
    hashmap_dna_count: &mut HashMap<usize, HashMap<DNA, usize>>,
    char: &DNA,
) {
    let inner_map = hashmap_dna_count
        .entry(pos)
        .or_insert_with(create_a_dna_to_count_hashmap);
    *inner_map.entry(char.clone()).or_insert(0) += 1;
}
/// Get hashmap value
fn get_hashmap_value(
    hashmap_dna_count: &HashMap<usize, HashMap<DNA, usize>>,
    pos: usize,
    key: DNA,
) -> usize {
    let tmp_hasmap = hashmap_dna_count.get(&pos).expect("should not happened");
    *tmp_hasmap.get(&key).expect("should not happend")
}

/// [[[ Write stats files ]]]
fn write_fq_stat_to_file(
    path: &str,
    hashmap_dna_count: &HashMap<usize, HashMap<DNA, usize>>,
    hashmap_low_qual_dna_count: &HashMap<usize, HashMap<DNA, usize>>,
    reads_size: &[u32],
    uniq_qual_value_count: usize,
) -> io::Result<()> {
    let col_name = "POS\t#bases\t%A\t%C\t%G\t%T\t%N\t%low\t%A\t%C\t%G\t%T";
    let mut output = File::create(path)?;
    writeln!(
        output,
        "[Quality value] {} distinct quality values.",
        uniq_qual_value_count
    )?;
    writeln!(output, "{}", get_read_len_stats_string(reads_size))?;
    writeln!(output, "{}", col_name)?;
    for pos in 0..(hashmap_dna_count.len()) {
        writeln!(
            output,
            "{}",
            get_dna_stats_string(hashmap_dna_count, hashmap_low_qual_dna_count, pos)
        )?;
    }
    Ok(())
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
    hashmap_dna_count: &HashMap<usize, HashMap<DNA, usize>>,
    hashmap_low_qual_dna_count: &HashMap<usize, HashMap<DNA, usize>>,
    pos: usize,
) -> String {
    let total: usize = hashmap_dna_count
        .get(&pos)
        .expect("should not happen")
        .values()
        .sum();
    let low_qual_total: usize = hashmap_low_qual_dna_count
        .get(&pos)
        .expect("should not happen")
        .values()
        .sum();
    let out = format!(
        "{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}",
        total,
        get_dna_proportion_string(hashmap_dna_count, pos, DNA::A, total),
        get_dna_proportion_string(hashmap_dna_count, pos, DNA::C, total),
        get_dna_proportion_string(hashmap_dna_count, pos, DNA::G, total),
        get_dna_proportion_string(hashmap_dna_count, pos, DNA::T, total),
        get_dna_proportion_string(hashmap_dna_count, pos, DNA::N, total),
        100.0 * low_qual_total as f64 / total as f64,
        get_dna_proportion_string(hashmap_low_qual_dna_count, pos, DNA::A, low_qual_total),
        get_dna_proportion_string(hashmap_low_qual_dna_count, pos, DNA::C, low_qual_total),
        get_dna_proportion_string(hashmap_low_qual_dna_count, pos, DNA::G, low_qual_total),
        get_dna_proportion_string(hashmap_low_qual_dna_count, pos, DNA::T, low_qual_total),
    );
    match pos {
        0 => {
            format!("Total\t{}", out)
        }
        _ => {
            format!("{}\t{}", pos, out)
        }
    }
}
fn get_dna_proportion_string(
    hashmap_dna_count: &HashMap<usize, HashMap<DNA, usize>>,
    pos: usize,
    dna: DNA,
    total: usize,
) -> f64 {
    100.0 * get_hashmap_value(hashmap_dna_count, pos, dna) as f64 / total as f64
}
