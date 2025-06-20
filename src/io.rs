use crate::bed::BedMap;

use bio::io::{fasta, fastq};

use flate2::read::GzDecoder;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter};

pub fn get_bed_map(file_path: &str) -> Result<BedMap, io::Error> {
    // Read the file and collect unique intervals for each name
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut bed_map: BedMap = BedMap::new();
    for line in reader.lines() {
        match line {
            Ok(line_content) => {
                let columns: Vec<&str> = line_content.split('\t').collect();
                if columns.len() > 2 {
                    let name = columns[0].to_string();
                    if let (Ok(start), Ok(end)) =
                        (columns[1].parse::<usize>(), columns[2].parse::<usize>())
                    {
                        bed_map.add(name, [start, end]);
                    } else {
                        eprintln!("Error parsing start or end for line: {}", line_content);
                    }
                }
            }
            Err(e) => eprintln!("Error reading line: {}", e),
        }
    }
    bed_map.merge();
    Ok(bed_map)
}
pub fn new_fq_iterator(fq_path: &str) -> io::Result<fastq::Reader<BufReader<Box<dyn BufRead>>>> {
    let reader: Box<dyn BufRead> = buffer_reader_maybe_gz(fq_path)?;
    Ok(fastq::Reader::new(reader))
}
pub fn new_fa_iterator(fa_path: &str) -> io::Result<fasta::Reader<BufReader<Box<dyn BufRead>>>> {
    let reader: Box<dyn BufRead> = buffer_reader_maybe_gz(fa_path)?;
    Ok(fasta::Reader::new(reader))
}
fn buffer_reader_maybe_gz(path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        let gz_decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(gz_decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}
pub fn append_bufwriter(out_path: &str) -> io::Result<BufWriter<File>> {
    let file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Open the file in append mode (no truncation)
        .open(out_path)?; // Open the file at the provided path
    let writer = BufWriter::new(file);
    Ok(writer)
}
