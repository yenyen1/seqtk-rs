use flate2::read::MultiGzDecoder;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, BufRead};
use std::path::Path;

enum FxType {
    Fasta,
    Fastq,
}
pub struct FxIterator {
    reader: Box<dyn BufRead>,
    file_type: FxType,
}

fn is_filetype(file_path: &Path, f_type: &str) -> bool {
    file_path.extension().and_then(OsStr::to_str) == Some(f_type)
}

fn get_fx_type(f_type: &str) -> Result<FxType, io::Error> {
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
    Ok(fx_type)
}

impl FxIterator {
    pub fn new(file_path: &str, f_type: &str) -> Result<FxIterator, io::Error> {
        let path = Path::new(file_path);
        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = if is_filetype(path, "gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        Ok(FxIterator {
            reader,
            file_type: get_fx_type(f_type).unwrap(),
        })
    }
}

unsafe impl Send for FxIterator {}

impl Iterator for FxIterator {
    type Item = Vec<String>;

    fn next(&mut self) -> Option<Self::Item> {
        let step = match self.file_type {
            FxType::Fasta => 2,
            FxType::Fastq => 4,
        };
        let mut fx_lines = vec![String::new(); step];
        for buf in fx_lines.iter_mut().take(step) {
            match self.reader.read_line(buf) {
                Ok(0) => {
                    return None;
                } // EOF reached, return None if no lines were read
                Ok(_) => {}
                Err(_) => {
                    return None;
                }
            }
        }
        Some(fx_lines)
    }
    
}
