use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

#[derive(Debug)]
enum DNA {
    A(bool),
    T(bool),
    C(bool),
    G(bool),
    N,
}

struct Read {
    name: String,
    seq: Vec<DNA>,
    qual: Vec<u8>,
}

impl Read {
    pub fn empty_read() -> Read {
        Read {
            name: String::new(),
            seq: Vec::new(),
            qual: Vec::new(),
        }
    }
    pub fn set_name(&mut self, name: &str) {
        self.name = name.split_whitespace().next().unwrap_or("").to_string();
    }
    pub fn set_seq(&mut self, seq_string: &str) {
        self.seq = seq_string
            .chars()
            .map(|c| match c {
                'A' => DNA::A(true),
                'a' => DNA::A(false),
                'T' => DNA::T(true),
                't' => DNA::T(false),
                'C' => DNA::C(true),
                'c' => DNA::C(false),
                'G' => DNA::G(true),
                'g' => DNA::G(false),
                _ => DNA::N,
            })
            .collect();
    }
    pub fn set_qual(&mut self, qual_string: &str, ascii_base: u8) {
        match ascii_base {
            33 | 64 => {
                self.qual = qual_string.chars().map(|c| c as u8 - ascii_base).collect();
            }
            _ => {
                println!("ascii base not support");
            }
        }
    }
    pub fn get_seq(&self) -> &Vec<DNA> {
        &self.seq
    }
    pub fn get_name(&self) -> &str {
        self.name.as_str()
    }
    pub fn get_qual(&self) -> &Vec<u8> {
        &self.qual
    }
    pub fn get_seq_len(&self) -> usize {
        self.seq.len()
    }
}

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

    let mut cur_read = Read::empty_read();
    let mut reads_size = Vec::<usize>::new();
    for (i, line) in reader.lines().enumerate() {
        let line_string = line.expect("input line is not a String");
        let first_char = line_string.get(0..1).expect("line_string out of index");
        match i % 4 {
            0 => {
                if first_char == "@" {
                    cur_read = Read::empty_read();
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
    println!("{}", get_read_stats(&mut reads_size));
    Ok(())
}

pub fn get_read_stats(reads_size: &mut [usize]) -> String {
    reads_size.sort();
    format!(
        "Length: mean={} ; min={} ; med={} ; max={} ; N50={}",
        average(reads_size),
        reads_size.iter().min().unwrap(),
        median(reads_size),
        reads_size.iter().max().unwrap(),
        n50(reads_size)
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

pub fn median(arr: &[usize]) -> f64 {
    let len_a = arr.len();
    match len_a {
        0 => 0.0,
        _ => {
            let mut tmp_arr = arr.to_vec();
            tmp_arr.sort();
            let mid = len_a / 2;
            match len_a % 2 {
                1 => tmp_arr[mid] as f64,
                _ => ((tmp_arr[mid - 1] as f64) + (tmp_arr[mid] as f64)) / 2.0,
            }
        }
    }
}

pub fn average(arr: &[usize]) -> f64 {
    let len_a = arr.len();
    match len_a {
        0 => 0.0,
        _ => arr.iter().sum::<usize>() as f64 / len_a as f64,
    }
}
pub fn n50(arr: &[usize]) -> usize {
    let len_a = arr.len();
    match len_a {
        0 => 0,
        _ => {
            let mut tmp_arr = arr.to_vec();
            tmp_arr.sort();
            let half_sum = tmp_arr.iter().sum::<usize>() / 2;
            let mut cur_acc_bases = 0;
            let mut cur_read_size = 0;
            while cur_acc_bases <= half_sum {
                cur_read_size = tmp_arr.pop().unwrap();
                cur_acc_bases += cur_read_size;
            }
            cur_read_size
        }
    }
}
