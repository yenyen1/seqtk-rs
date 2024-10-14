use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

enum FxType {
    Fasta,
    Fastq,
}
pub struct FxIterator {
    reader: Box<dyn BufRead>,
    file_type: FxType,
}

impl FxIterator {
    pub fn new(file_path: &str, f_type: &str, is_gz: bool) -> Result<FxIterator, io::Error> {
        let file = File::open(file_path)?;
        let reader: Box<dyn BufRead> = if is_gz {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        let fx_type = match f_type {
            "fastq" | "fq" => FxType::Fastq,
            "fasta" | "fa" => FxType::Fasta,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "wrong file extension.",
                ))
            }
        };

        Ok(FxIterator {
            reader: reader,
            file_type: fx_type,
        })
    }
}

impl Iterator for FxIterator {
    type Item = Vec<String>;

    fn next(&mut self) -> Option<Self::Item> {
        let step = match self.file_type {
            FxType::Fasta => 2,
            FxType::Fastq => 4,
        };
        let mut fx_lines = vec![String::new(); step];
        for i in 0..step {
            match self.reader.read_line(&mut fx_lines[i]) {
                Ok(_) => {}
                Err(_) => {
                    return None;
                }
            }
        }
        Some(fx_lines)
    }
}
