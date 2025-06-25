use crate::dna::SeqComp;
use crate::io_utils::{FaReader, FqReader, Output};
use std::io;

pub fn calc_fq_comp(path: &str) -> Result<(), std::io::Error> {
    let fq_iter = FqReader::new(path)?;
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    for record in fq_iter.records() {
        let read = record.unwrap();
        let count = SeqComp::count_nucleotides(read.seq(), 0, read.seq().len());
        output.write(format!(
            "{}\n",
            count
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        ))?;
    }
    Ok(())
}
pub fn calc_fa_comp(path: &str) -> Result<(), std::io::Error> {
    let fa_iter = FaReader::new(path)?;
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    for record in fa_iter.records() {
        let read = record.unwrap();
        let count = SeqComp::count_nucleotides(read.seq(), 0, read.seq().len());
        output.write(format!(
            "{}\n",
            count
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        ))?;
    }
    Ok(())
}
