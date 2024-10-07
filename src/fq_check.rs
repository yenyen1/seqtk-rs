use crate::sequence_read::SequenceRead;
use crate::stats;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

pub fn fq_check(fq_path: &str, quality_value: Option<u8>) {
    let fq_post_fix = fq_path.split(".").last().unwrap_or("");
    match fq_post_fix {
        "fq" | "fastq" => {
            process_fastq_and_check(fq_path, quality_value)
                .expect("[Error] failed to process fastq file.");
        }
        "fa" | "fasta" => {
            println!("add fasta later")
        }
        "gz" => {
            let fq_post_fix_wo_gz = fq_path.rsplit(".").nth(1).unwrap_or("");
            match fq_post_fix_wo_gz {
                "fq" | "fastq" => {
                    process_fastq_gz(fq_path, quality_value)
                        .expect("[Error] failed to process fastq.gz file.");
                }
                "fa" | "fasta" => {
                    println!("add fasta later")
                }
                _ => {
                    println!("[Error] Input file format {} not supported", fq_post_fix);
                }
            }
        }
        _ => {
            println!("[Error] Input file format {} not supported", fq_post_fix);
        }
    }
}

pub fn process_fastq_and_check(fq_path: &str, quality_value: Option<u8>) -> io::Result<()> {
    // let q = quality_value.unwrap_or(0);
    // println!("process fastq: {} {}", fq_path, q);

    let file = File::open(fq_path)?;
    let reader = BufReader::new(file);

    let mut cur_read = SequenceRead::empty_read();
    let mut reads_size = Vec::<usize>::new();
    for (i, line) in reader.lines().enumerate() {
        let line_string = line.expect("input line is not a String");
        let first_char = line_string.get(0..1).expect("line_string out of index");
        match i % 4 {
            0 => {
                if first_char == "@" {
                    cur_read = SequenceRead::empty_read();
                    cur_read.set_name(line_string.get(1..).expect("line_string out of index"));
                } else {
                    println!("Not found read name char @ in fastq file");
                }
            }
            1 => match first_char {
                "A" | "T" | "C" | "G" | "N" => {
                    cur_read.set_seq(&line_string);
                }
                _ => {
                    println!("Not an expected char (ATCGN) for the first char in seq.");
                }
            },
            3 => {
                cur_read.set_qual(&line_string, 33);
                reads_size.push(cur_read.get_seq_len());
                // println!("test name: {}", cur_read.get_name());
                // println!("test seq len: {}", cur_read.get_seq_len());
                // println!("test qual: {:?}", cur_read.get_qual());
                // break;
            }
            _ => {}
        }
    }
    println!("{}", get_read_stats(&reads_size));
    Ok(())
}

pub fn get_read_stats(reads_size: &[usize]) -> String {
    let mut sorted_reads_size = reads_size.to_vec();
    sorted_reads_size.sort();
    format!(
        "Length: mean={} ; min={} ; med={} ; max={} ; N50={}",
        stats::average(&sorted_reads_size),
        sorted_reads_size.iter().min().unwrap(),
        stats::median(&sorted_reads_size),
        sorted_reads_size.iter().max().unwrap(),
        stats::n50(&sorted_reads_size)
    )
}

pub fn process_fastq_gz(fq_path: &str, quality_value: Option<u8>) -> io::Result<()> {
    let q = quality_value.unwrap_or(0);
    println!("process fastq: {} {}", fq_path, q);

    let file = File::open(fq_path)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    for line in reader.lines() {
        let line_string = line.expect("not a line");
        println!("print > {}", line_string);
    }
    Ok(())
}
