use crate::bed::BedMap;
use crate::dna::SeqComp;
use crate::io_utils::{FaReader, FqReader, Output};
use std::io::{self, Stdout};

pub fn calc_fq_comp_wo_bed(path: &str, exclude_masked: bool) -> Result<(), std::io::Error> {
    let fq_iter = FqReader::new(path)?;
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    print_cols(&mut output)?;
    if exclude_masked {
        for record in fq_iter.records() {
            let read = record.unwrap();
            let mut count: [usize; 23] = [0; 23];
            SeqComp::count_unmasked_nc(&mut count, read.seq(), 0, read.seq().len());
            let result = SeqComp::get_unmasked_result(&count);
            print(&mut output, read.id(), read.seq().len(), &result)?;
        }
    } else {
        for record in fq_iter.records() {
            let read = record.unwrap();
            let mut count: [usize; 23] = [0; 23];
            SeqComp::count_all_nc(&mut count, read.seq(), 0, read.seq().len());
            let result = SeqComp::get_all_result(&count);
            print(&mut output, read.id(), read.seq().len(), &result)?;
        }
    }
    Ok(())
}
pub fn calc_fa_comp_wo_bed(path: &str, exclude_masked: bool) -> Result<(), std::io::Error> {
    let fa_iter = FaReader::new(path)?;
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    print_cols(&mut output)?;
    if exclude_masked {
        for record in fa_iter.records() {
            let read = record.unwrap();
            let mut count: [usize; 23] = [0; 23];
            SeqComp::count_unmasked_nc(&mut count, read.seq(), 0, read.seq().len());
            let result = SeqComp::get_unmasked_result(&count);
            print(&mut output, read.id(), read.seq().len(), &result)?;
        }
    } else {
        for record in fa_iter.records() {
            let read = record.unwrap();
            let mut count: [usize; 23] = [0; 23];
            SeqComp::count_all_nc(&mut count, read.seq(), 0, read.seq().len());
            let result = SeqComp::get_all_result(&count);
            print(&mut output, read.id(), read.seq().len(), &result)?;
        }
    }
    Ok(())
}

pub fn calc_fq_comp_with_bed(
    path: &str,
    bed: &str,
    exclude_masked: bool,
) -> Result<(), std::io::Error> {
    let fq_iter = FqReader::new(path)?;
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    print_cols(&mut output)?;
    if exclude_masked {
        let bedmap = BedMap::from(bed)?;
        for record in fq_iter.records() {
            let read = record.unwrap();
            let mut count: [usize; 23] = [0; 23];
            if let Some(bedvec) = bedmap.get(read.id()) {
                bedvec.iter().for_each(|pos| {
                    SeqComp::count_unmasked_nc(&mut count, read.seq(), pos.0, pos.1);
                });
                let result = SeqComp::get_unmasked_result(&count);
                print(&mut output, read.id(), read.seq().len(), &result)?;
            }
        }
    } else {
        let bedmap = BedMap::from(bed)?;
        for record in fq_iter.records() {
            let read = record.unwrap();
            let mut count: [usize; 23] = [0; 23];
            if let Some(bedvec) = bedmap.get(read.id()) {
                bedvec.iter().for_each(|pos| {
                    SeqComp::count_all_nc(&mut count, read.seq(), pos.0, pos.1);
                });
                let result = SeqComp::get_all_result(&count);
                print(&mut output, read.id(), read.seq().len(), &result)?;
            }
        }
    }
    Ok(())
}

pub fn calc_fa_comp_with_bed(
    path: &str,
    bed: &str,
    exclude_masked: bool,
) -> Result<(), std::io::Error> {
    let fa_iter = FaReader::new(path)?;
    let stdout = io::stdout();
    let mut output = Output::new(stdout);
    print_cols(&mut output)?;
    if exclude_masked {
        let bedmap = BedMap::from(bed)?;
        for record in fa_iter.records() {
            let read = record.unwrap();
            let mut count: [usize; 23] = [0; 23];
            if let Some(bedvec) = bedmap.get(read.id()) {
                bedvec.iter().for_each(|pos| {
                    SeqComp::count_unmasked_nc(&mut count, read.seq(), pos.0, pos.1);
                });
                let result = SeqComp::get_unmasked_result(&count);
                print(&mut output, read.id(), read.seq().len(), &result)?;
            }
        }
    } else {
        let bedmap = BedMap::from(bed)?;
        for record in fa_iter.records() {
            let read = record.unwrap();
            let mut count: [usize; 23] = [0; 23];
            if let Some(bedvec) = bedmap.get(read.id()) {
                bedvec.iter().for_each(|pos| {
                    SeqComp::count_all_nc(&mut count, read.seq(), pos.0, pos.1);
                });
                let result = SeqComp::get_all_result(&count);
                print(&mut output, read.id(), read.seq().len(), &result)?;
            }
        }
    }
    Ok(())
}

fn print_cols(output: &mut Output<Stdout>) -> Result<(), std::io::Error> {
    output.write("ID\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CG\t#GC\n")?;
    Ok(())
}
fn print(
    output: &mut Output<Stdout>,
    id: &str,
    size: usize,
    count: &[usize; 9],
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
