use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

pub fn fq_check(fq_path: &str, quality_value: Option<u32>) {
    let fq_path_split_vec = fq_path.rsplit(".").collect::<Vec<_>>();
    let fq_post_fix = fq_path_split_vec.first().expect("out of index");
    match *fq_post_fix {
        "fq" | "fastq" => {
            let _ = process_fastq(fq_path, quality_value);
        }
        "gz" => {
            let _ = process_fastq_gz(fq_path, quality_value);
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
    for line in reader.lines() {
        let line_string = line.expect("not a line");
        println!("print > {}", line_string);
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
