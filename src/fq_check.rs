use crate::dna;
use crate::sequence_read::SequenceRead;
use crate::stats;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
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
    let mut hashmap_db: HashMap<usize, HashMap<char, usize>> =
        HashMap::from([(0, init_dna_char_hashmap())]);
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
                add_sequence_count_into_hashmap(&mut hashmap_db, &cur_read);
            }
            _ => {}
        }
    }
    // println!("{}", get_read_stats_string(&reads_size));
    write_fq_stat_to_file("test.txt", &hashmap_db, &reads_size)?;

    Ok(())
}

pub fn write_fq_stat_to_file(
    path: &str,
    hashmap_db: &HashMap<usize, HashMap<char, usize>>,
    reads_size: &[usize],
) -> io::Result<()> {
    let mut output = File::create(path)?;
    write!(output, "{}", get_read_stats_string(reads_size))?;
    // for (key, val) in hashmap_db.iter() {
    //     println!("{} : {:?}", key, val);
    // }
    Ok(())
}

pub fn add_sequence_count_into_hashmap(
    hashmap_db: &mut HashMap<usize, HashMap<char, usize>>,
    read: &SequenceRead,
) {
    for (i, char) in read.get_seq().iter().enumerate() {
        add_dna_char_count_into_pos(0, hashmap_db, char);
        add_dna_char_count_into_pos(i + 1, hashmap_db, char);
    }
}

pub fn add_dna_char_count_into_pos(
    pos: usize,
    hashmap_db: &mut HashMap<usize, HashMap<char, usize>>,
    char: &dna::DNA,
) {
    match char {
        dna::DNA::A(true) | dna::DNA::A(false) => {
            add_char_count_into_pos(pos, hashmap_db, 'A');
        }
        dna::DNA::T(true) | dna::DNA::T(false) => {
            add_char_count_into_pos(pos, hashmap_db, 'T');
        }
        dna::DNA::C(true) | dna::DNA::C(false) => {
            add_char_count_into_pos(pos, hashmap_db, 'C');
        }
        dna::DNA::G(true) | dna::DNA::G(false) => {
            add_char_count_into_pos(pos, hashmap_db, 'G');
        }
        dna::DNA::N => {
            add_char_count_into_pos(pos, hashmap_db, 'N');
        }
    }
}

pub fn add_char_count_into_pos(
    pos: usize,
    hashmap_db: &mut HashMap<usize, HashMap<char, usize>>,
    char: char,
) {
    hashmap_db
        .entry(pos)
        .or_insert(init_dna_char_hashmap())
        .entry(char)
        .and_modify(|c| *c += 1)
        .or_insert(1);
}

pub fn init_dna_char_hashmap() -> HashMap<char, usize> {
    HashMap::from([('A', 0), ('T', 0), ('C', 0), ('G', 0), ('N', 0)])
}

pub fn get_read_stats_string(reads_size: &[usize]) -> String {
    let mut sorted_reads_size = reads_size.to_vec();
    sorted_reads_size.sort();
    format!(
        "Length: mean={} ; min={} ; med={} ; max={} ; N50={}\n",
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
