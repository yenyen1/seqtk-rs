use crate::io::recordset::OwnedRecordSet;

use std::fs::File;
use std::io::{Error, ErrorKind, BufRead, BufReader, Read};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;
use seq_io::fasta::{self, Record as _};
use seq_io::fastq::{self, Record};

type BoxedRead = Box<dyn Read + Send>;

pub enum FxReader {
    Fasta(fasta::Reader<BoxedRead>),
    Fastq(fastq::Reader<BoxedRead>),
}
impl FxReader {
    pub fn new<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        let is_gzip = {
            let buf = reader.fill_buf()?;
            buf.starts_with(&[0x1f, 0x8b])
        };

        let decoder: BoxedRead = if is_gzip {
            Box::new(MultiGzDecoder::new(reader))
        } else {
            Box::new(reader)
        };

        let mut format_reader = BufReader::new(decoder);
        let first_byte = {
            let buf = format_reader.fill_buf()?;
            if buf.is_empty() {return Err(Error::new(ErrorKind::UnexpectedEof, "Empty"));}
            buf[0]
        };

        let boxed_reader: BoxedRead = Box::new(format_reader);
        match first_byte {
            b'>' => Ok(FxReader::Fasta(fasta::Reader::new(boxed_reader))),
            b'@' => Ok(FxReader::Fastq(fastq::Reader::new(boxed_reader))),
            _ => Err(Error::new(ErrorKind::InvalidData, "Unknown format")),
        }
    } 
}

pub struct BatchReader {
    inner: FxReader,
    leftover_fa: Option<fasta::OwnedRecord>,
    leftover_fq: Option<fastq::OwnedRecord>,
}
impl BatchReader {
    pub fn new(inner: FxReader) -> Self {
        Self {
            inner,
            leftover_fa: None,
            leftover_fq: None,
        }
    }

    pub fn fill_batch(&mut self, batch: &mut OwnedRecordSet) -> std::io::Result<bool> {
        batch.clear();

        match &mut self.inner {
            FxReader::Fasta(reader) => {
                if let Some(record) = self.leftover_fa.take() {
                    batch.push(record.id_bytes(), record.seq(), &[]);
                }

                while let Some(record) = reader.next() {
                    match record {
                        Ok(rec) => {
                            let seq = rec.full_seq();
                            if batch.is_overload(seq.len()) {
                                self.leftover_fa = Some(rec.to_owned_record());
                                break;
                            }
                            batch.push(rec.id_bytes(), &seq, &[]);
                        },
                        Err(e) => { 
                            eprintln!("FASTA parsing Error (Skip record): {}", e);
                        },
                    };
                };
            },

            FxReader::Fastq(reader) => {
                if let Some(record) = self.leftover_fq.take() {
                    batch.push(record.id_bytes(), record.seq(), &[]);
                }

                while let Some(record) = reader.next() {
                    match record {
                        Ok(rec) => {
                            if batch.is_overload(rec.seq().len()) {
                                self.leftover_fq = Some(rec.to_owned_record());
                            }
                            batch.push(rec.id_bytes(), rec.seq(), rec.qual());
                        },
                        Err(e) => {
                            eprintln!("FASTQ parsing Error (Skip record): {}", e);
                        },
                    }
                }
            },
        }
        Ok(!batch.is_empty()) // return True if batch read something
    }
}