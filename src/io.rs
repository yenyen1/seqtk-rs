
use bio::io::{fasta, fastq};
use flate2::read::GzDecoder;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Stdout};

// pub fn new_fq_iterator(fq_path: &str) -> io::Result<fastq::Reader<BufReader<Box<dyn BufRead>>>> {
//     let reader: Box<dyn BufRead> = buffer_reader_maybe_gz(fq_path)?;
//     Ok(fastq::Reader::new(reader))
// }
// pub fn new_fa_iterator(fa_path: &str) -> io::Result<fasta::Reader<BufReader<Box<dyn BufRead>>>> {
//     let reader: Box<dyn BufRead> = buffer_reader_maybe_gz(fa_path)?;
//     Ok(fasta::Reader::new(reader))
// }
pub fn buffer_reader_maybe_gz(path: &str) -> io::Result<Box<dyn BufRead>> {
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

pub struct FaReader(fasta::Reader<BufReader<Box<dyn BufRead>>>);
pub struct FqReader(fastq::Reader<BufReader<Box<dyn BufRead>>>);
impl FaReader {
    pub fn new(path: &str) -> io::Result<Self> {
        let reader: Box<dyn BufRead> = buffer_reader_maybe_gz(path)?;
        Ok(Self(fasta::Reader::new(reader)))
    }
    pub fn records(self) -> fasta::Records<BufReader<Box<dyn BufRead + 'static>>> {
        self.0.records()
    }
}
impl FqReader {
    pub fn new(path: &str) -> io::Result<Self> {
        let reader: Box<dyn BufRead> = buffer_reader_maybe_gz(path)?;
        Ok(Self(fastq::Reader::new(reader)))
    }
    pub fn records(self) -> fastq::Records<BufReader<Box<dyn BufRead + 'static>>> {
        self.0.records()
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
