use bio::io::fastq;
use flate2::read::GzDecoder;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter};
use std::path::Path;

pub fn append_bufwriter(out_path: &str) -> io::Result<BufWriter<File>> {
    let file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Open the file in append mode (no truncation)
        .open(out_path)?; // Open the file at the provided path
    let writer = BufWriter::new(file);
    Ok(writer)
}
pub fn new_fx_iterator(file_path: &str) -> io::Result<fastq::Reader<BufReader<Box<dyn BufRead>>>> {
    let file_extension = Path::new(file_path).extension().and_then(|s| s.to_str());
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = match file_extension {
        Some("gz") => {
            let gz_decoder = GzDecoder::new(file);
            Box::new(BufReader::new(gz_decoder))
        }
        _ => Box::new(BufReader::new(file)),
    };
    Ok(fastq::Reader::new(reader))
}
