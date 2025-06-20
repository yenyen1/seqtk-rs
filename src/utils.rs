use bio::io::{fasta, fastq};
use flate2::read::GzDecoder;

use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Stdout};
use std::path::Path;

pub fn append_bufwriter(out_path: &str) -> io::Result<BufWriter<File>> {
    let file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Open the file in append mode (no truncation)
        .open(out_path)?; // Open the file at the provided path
    let writer = BufWriter::new(file);
    Ok(writer)
}
pub fn new_fq_iterator(file_path: &str) -> io::Result<fastq::Reader<BufReader<Box<dyn BufRead>>>> {
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
pub fn new_fa_iterator(file_path: &str) -> io::Result<fasta::Reader<BufReader<Box<dyn BufRead>>>> {
    let file_extension = Path::new(file_path).extension().and_then(|s| s.to_str());
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = match file_extension {
        Some("gz") => {
            let gz_decoder = GzDecoder::new(file);
            Box::new(BufReader::new(gz_decoder))
        }
        _ => Box::new(BufReader::new(file)),
    };
    Ok(fasta::Reader::new(reader))
}
pub trait RecordType {
    fn seq(&self) -> &[u8];
    fn id(&self) -> &str;
    fn desc(&self) -> Option<&str>;
    fn qual(&self) -> &[u8];
}
impl RecordType for fasta::Record {
    fn seq(&self) -> &[u8] {
        self.seq()
    }
    fn id(&self) -> &str {
        self.id()
    }
    fn desc(&self) -> Option<&str> {
        self.desc()
    }
    fn qual(&self) -> &[u8] {
        panic!("[RecordType::Fasta] no qual function.")
    }
}
impl RecordType for fastq::Record {
    fn seq(&self) -> &[u8] {
        self.seq()
    }
    fn id(&self) -> &str {
        self.id()
    }
    fn desc(&self) -> Option<&str> {
        self.desc()
    }
    fn qual(&self) -> &[u8] {
        self.qual()
    }
}
pub enum FxWriter {
    Fasta(fasta::Writer<Stdout>),
    Fastq(fastq::Writer<Stdout>),
}
impl FxWriter {
    pub fn new(is_fasta: bool) -> FxWriter {
        let std_out = io::stdout();
        if is_fasta {
            FxWriter::Fasta(fasta::Writer::new(std_out))
        } else {
            FxWriter::Fastq(fastq::Writer::new(std_out))
        }
    }
    pub fn write(
        &mut self,
        id: &str,
        seq: &[u8],
        desc: Option<&str>,
        qual: &[u8],
    ) -> io::Result<()> {
        match self {
            FxWriter::Fasta(ref mut out) => {
                out.write(id, desc, seq)?;
            }
            FxWriter::Fastq(ref mut out) => {
                out.write(id, desc, seq, qual)?;
            }
        }
        Ok(())
    }
}
