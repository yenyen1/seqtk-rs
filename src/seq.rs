// use std::io;
use crate::utils;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

use bio::io::fastq::{self, Record};
// use std::io::{self, BufWriter, Write};

pub struct FilterRule {
    mini_seq_length: usize,
    drop_ambigous_seq: bool,
    output_odd_reads: bool,
    output_even_reads: bool,
    random_seed: u64,
    sample_fraction: Option<f64>,
}
impl FilterRule {
    pub fn new(
        mini_seq_length: usize,
        drop_ambigous_seq: bool,
        output_odd_reads: bool,
        output_even_reads: bool,
        random_seed: u64,
        sample_fraction: Option<f64>,
    ) -> Self {
        assert!(
            !(output_even_reads && output_odd_reads),
            "Should not use --output-even-reads and --output-odd-reads option together."
        );
        FilterRule {
            mini_seq_length,
            drop_ambigous_seq,
            output_odd_reads,
            output_even_reads,
            random_seed,
            sample_fraction,
        }
    }
}

pub struct ParsedSeqArgs {
    q_low: u8,
    q_high: u8,
    mask_char: Option<char>,
    mask_regions: Option<String>,
    quality_shift: u8,
    ///ascii_bases
    shift_quality_33: bool,
    uppercases: bool,
    lowercases: bool,

    n_residues: Option<u32>,

    fake_fastq_quality: bool,
    mask_complement_region: bool,
    reverse_complement: bool,
    both_complement: bool,
    output_fasta: bool,
    trim_header: bool,

    strip_whitespace: bool,
}

pub fn parse_seq(
    fx_path: &String,
    out_path: &String,
    filter_rule: &FilterRule,
) -> Result<(), std::io::Error> {
    match filter_rule.sample_fraction {
        Some(frac) => {
            let mut rng = StdRng::seed_from_u64(filter_rule.random_seed);
            let fx_iter = utils::new_fx_iterator(fx_path)?;
            for (i, record) in fx_iter.records().enumerate() {
                let read = record.unwrap();
                let filter = filter_read(i + 1, &read, filter_rule);
                if filter && rng.gen::<f64>() <= frac {
                    println!(
                        "[i={}] read length {} > {}",
                        i + 1,
                        read.seq().len(),
                        filter_rule.mini_seq_length
                    );
                    // break
                }
            }
        }
        None => {
            let fx_iter = utils::new_fx_iterator(fx_path)?;
            for (i, record) in fx_iter.records().enumerate() {
                let read = record.unwrap();
                let filtered_read = filter_read(i + 1, &read, filter_rule);
                if filtered_read {
                    println!(
                        "[i={}] read length {} > {}",
                        i + 1,
                        read.seq().len(),
                        filter_rule.mini_seq_length
                    );
                    // break
                }
            }
        }
    }
    Ok(())
}

fn modify_record(read: &Record) {
    //add
}

fn filter_read(i: usize, read: &Record, filter_rule: &FilterRule) -> bool {
    let even_or_odd = if !(filter_rule.output_even_reads && filter_rule.output_odd_reads) {
        if filter_rule.output_even_reads {
            i % 2 == 0
        } else if filter_rule.output_odd_reads {
            i % 2 == 1
        } else {
            true
        }
    } else {
        false
    };
    let seq_wo_n: bool = if filter_rule.drop_ambigous_seq {
        !read.seq().contains(&078_u8) // seq not containing N
    } else {
        true
    };
    let write_read = read.seq().len().gt(&filter_rule.mini_seq_length) && // seq len large than min_len 
        seq_wo_n && // seq not containing N
        even_or_odd;
    write_read
}
