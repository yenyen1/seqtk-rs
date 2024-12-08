// use std::io;
use crate::utils;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

use bio::io::fastq::{self, Record};
// use std::io::{self, BufWriter, Write};

pub struct FilterParas {
    mini_seq_length: usize,
    drop_ambigous_seq: bool,
    output_odd_reads: bool,
    output_even_reads: bool,
    random_seed: u64,
    sample_fraction: Option<f64>,
}
impl FilterParas {
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
            "Should not use --output-even-reads and --output-odd-reads options together."
        );
        FilterParas {
            mini_seq_length,
            drop_ambigous_seq,
            output_odd_reads,
            output_even_reads,
            random_seed,
            sample_fraction,
        }
    }
}

pub struct MaskParas<'a> {
    q_low: u8,
    q_high: u8,
    mask_char: Option<char>,
    mask_regions: &'a Option<String>,
    mask_complement_region: bool,
    uppercases: bool,
    lowercases_to_char: bool,
}
impl<'a> MaskParas<'a> {
    pub fn new(
        q_low: u8,
        q_high: u8,
        mask_char: Option<char>,
        mask_regions: &'a Option<String>,
        mask_complement_region: bool,
        uppercases: bool,
        lowercases_to_char: bool,
    ) -> Self {
        assert!(
            !(mask_complement_region && mask_regions.is_none()),
            "--mask-complment-region should effective with --mask-regions."
        );
        assert!(
            !(uppercases && (lowercases_to_char || mask_char.is_some())),
            "Should not use --uppercases with --lowercases_to_char or mask-char options together."
        );
        assert!(
            !(lowercases_to_char && mask_char.is_none()),
            "--lowercases-to-char should effective with --mask-char."
        );
        MaskParas {
            q_low,
            q_high,
            mask_char,
            mask_regions,
            mask_complement_region,
            uppercases,
            lowercases_to_char,
        }
    }
}
pub struct OutputArgs {
    quality_shift: u8, // ascii_bases
    shift_quality_33: bool,

    n_residues: Option<u32>,
    reverse_complement: bool,
    both_complement: bool,
    output_fasta: bool,
    fake_fastq_quality: bool,
    trim_header: bool,
    strip_whitespace: bool,
}

pub fn parse_seq(
    fx_path: &String,
    out_path: &String,
    filter_rule: &FilterParas,
    mask_paras: &MaskParas,
) -> Result<(), std::io::Error> {
    match filter_rule.sample_fraction {
        Some(frac) => {
            let mut rng = StdRng::seed_from_u64(filter_rule.random_seed);
            let fx_iter = utils::new_fx_iterator(fx_path)?;
            for (i, record) in fx_iter.records().enumerate() {
                let read = record.unwrap();
                let filter = filter_read(i + 1, &read, filter_rule);
                if filter && rng.gen::<f64>() <= frac {
                    let seq = modify_seq(&read, &mask_paras);
                    println!("[i={}] {:?}", i, String::from_utf8(seq));
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
                    let seq = modify_seq(&read, &mask_paras);
                    println!("[i={}] {:?}", i, String::from_utf8(seq));
                    // break
                }
            }
        }
    }
    Ok(())
}

fn modify_seq(read: &Record, mask_paras: &MaskParas) -> Vec<u8> {
    let mut seq = if mask_paras.uppercases {
        read.seq().to_ascii_uppercase()
    } else {
        read.seq().to_vec()
    };
    match mask_paras.mask_char {
        Some(c) => {
            // masked bases to char.
            let c_u8 = c as u8;
            if mask_paras.lowercases_to_char {
                for (&qual, ch) in read.qual().iter().zip(seq.iter_mut()) {
                    if qual < mask_paras.q_low || mask_paras.q_high < qual {
                        *ch = c_u8; // mask bases by q_low and q_high
                    } else if ch.is_ascii_lowercase() {
                        *ch = c_u8; // convert lowercase to char
                    }
                }
            } else {
                for (&qual, ch) in read.qual().iter().zip(seq.iter_mut()) {
                    if qual < mask_paras.q_low || mask_paras.q_high < qual {
                        *ch = c_u8; // mask bases by q_low and q_high
                    }
                }
            }
        }
        None => {
            // masked bases to lowercases.
            for (&qual, ch) in read.qual().iter().zip(seq.iter_mut()) {
                if qual < mask_paras.q_low || mask_paras.q_high < qual {
                    *ch = ch.to_ascii_lowercase(); // mask bases by q_low and q_high
                }
            }
        }
    }
    seq
}

fn filter_read(i: usize, read: &Record, filter_rule: &FilterParas) -> bool {
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
