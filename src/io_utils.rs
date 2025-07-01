use bio::io::{fasta, fastq};
use flate2::read::GzDecoder;
use std::fmt::Display;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Stdout, Write};

pub fn buffer_reader_maybe_gz(path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        let gz_decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(gz_decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
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

pub struct Output<W: Write> {
    writer: W,
}
impl<W: Write> Output<W> {
    pub fn new(writer: W) -> Self {
        Output { writer }
    }
    pub fn write<T: Display>(&mut self, result: T) -> io::Result<()> {
        write!(self.writer, "{}", result)
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
