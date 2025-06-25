use crate::dna::SeqComp;
use crate::io_utils::{FaReader, FqReader, Output};
use std::io::{self, Stdout};

pub fn calc_fq_comp(path: &str, exclude_masked: bool) -> Result<(), std::io::Error> {
    let fq_iter = FqReader::new(path)?;
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    print_cols(&mut output)?;
    if exclude_masked {
        for record in fq_iter.records() {
            let read = record.unwrap();
            let count = SeqComp::count_unmasked_nucleotides(read.seq(), 0, read.seq().len());
            print(&mut output, read.id(), read.seq().len(), &count)?;
        }
    } else {
        for record in fq_iter.records() {
            let read = record.unwrap();
            let count = SeqComp::count_all_nucleotides(read.seq(), 0, read.seq().len());
            print(&mut output, read.id(), read.seq().len(), &count)?;
        }
    }
    Ok(())
}
pub fn calc_fa_comp(path: &str, exclude_masked: bool) -> Result<(), std::io::Error> {
    let fa_iter = FaReader::new(path)?;
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    print_cols(&mut output)?;
    if exclude_masked {
        for record in fa_iter.records() {
            let read = record.unwrap();
            let count = SeqComp::count_unmasked_nucleotides(read.seq(), 0, read.seq().len());
            print(&mut output, read.id(), read.seq().len(), &count)?;
        }
    } else {
        for record in fa_iter.records() {
            let read = record.unwrap();
            let count = SeqComp::count_all_nucleotides(read.seq(), 0, read.seq().len());
            print(&mut output, read.id(), read.seq().len(), &count)?;
        }
    }
    Ok(())
}
fn print_cols(output: &mut Output<Stdout>) -> Result<(), std::io::Error> {
    output.write("ID\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts\n")?;
    Ok(())
}
fn print(
    output: &mut Output<Stdout>,
    id: &str,
    size: usize,
    count: &[usize; 11],
) -> Result<(), std::io::Error> {
    output.write(format!(
        "{}\t{}\t{}\n",
        id,
        size,
        count
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join("\t")
    ))?;
    Ok(())
}
