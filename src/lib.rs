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
}

pub fn fq_check(fq_path: &str, quality_value: Option<u32>) {
    let fq_path_split_vec = fq_path.rsplit(".").collect::<Vec<_>>();
    let fq_post_fix = fq_path_split_vec.first().expect("out of index");
    match *fq_post_fix {
        "fq" | "fastq" => {
            process_fastq(fq_path, quality_value).expect("[Error] failed to process fastq file.");
        }
        "fa" | "fasta" => {
            println!("add fasta later")
        }
        "gz" => {
            let fq_post_fix_wo_gz = fq_path_split_vec.get(1).expect("out of index");
            match *fq_post_fix_wo_gz {
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

pub fn process_fastq(fq_path: &str, quality_value: Option<u32>) -> io::Result<()> {
    let q = quality_value.unwrap_or(0);
    println!("process fastq: {} {}", fq_path, q);

    let file = File::open(fq_path)?;
    let reader = BufReader::new(file);
    let mut cur_read = Read::empty_read();
    for (i, line) in reader.lines().enumerate() {
        let line_string = line.expect("not a line");
        println!("print > {}", line_string);
        let first_char = line_string.get(0..1).expect("out of index");
        match i % 4 {
            0 => {
                if first_char == "@" {
                    cur_read = Read::empty_read();
                    cur_read.set_name(line_string.get(1..).expect("out of index"));
                } else {
                    println!("Not read name char @");
                }
            }
            1 => match first_char {
                "A" | "T" | "C" | "G" | "N" => {
                    cur_read.set_seq(&line_string);
                }
                _ => {
                    println!("Not expect seq char.");
                }
            },
            3 => {
                cur_read.set_qual(&line_string, 33);
                println!("test name: {}", cur_read.get_name());
                println!("test seq: {:?}", cur_read.get_seq());
                println!("test qual: {:?}", cur_read.get_qual());
                break;
            }
            _ => {}
        }
    }

    Ok(())
}

pub fn process_fastq_gz(fq_path: &str, quality_value: Option<u32>) -> io::Result<()> {
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
