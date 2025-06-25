use crate::io_utils::{FaReader, FqReader, Output};
use std::io;

pub fn calc_fq_size(path: &str) -> Result<(), std::io::Error> {
    let fq_iter = FqReader::new(path)?;
    let mut seq_count = 0;
    let mut bases_count = 0;
    for record in fq_iter.records() {
        let read = record.unwrap();
        seq_count += 1;
        bases_count += read.seq().len();
    }
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    output.write(format!("{}\t{}\n", seq_count, bases_count))?;
    Ok(())
}

pub fn calc_fa_size(path: &str) -> Result<(), std::io::Error> {
    let fa_iter = FaReader::new(path)?;
    let mut seq_count = 0;
    let mut bases_count = 0;
    for record in fa_iter.records() {
        let read = record.unwrap();
        seq_count += 1;
        bases_count += read.seq().len();
    }
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    output.write(format!("{}\t{}\n", seq_count, bases_count))?;
    Ok(())
}
