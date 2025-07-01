use crate::bed::BedMap;
use crate::dna::SeqComp;
use crate::io_utils::{FaReader, FqReader, Output};
use crate::record::RecordType;

/// Parses FASTQ file and compute the statistic w/o masked sequences.
/// Outputs the results to [`std::io::stdout()`].
///The output columns are:
/// - `#seq`: Number of reads
/// - `#bases`: Number of bases
/// - `#A`, `#C`, `#G`, `#T`: Number of each nucleotide
/// - `#2`: Number of `R`, `Y`, `S`, `W`, `K`, `M`
/// - `#3`: Number of `B`, `D`, `H`, `V`
/// - `#4`: Number of `N`
/// - `CG`: Number of `CG` on the template strand
/// - `GC`: Number of `GC` on the template strand
///
/// # Arguments
///
/// * `path` - FASTQ path
/// * `exclude_masked` - If true, masked sequences (e.g., lowercases or char other then IUPAC code) will be excluded from the output.
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn calc_fq_comp_wo_bed(path: &str, exclude_masked: bool) -> Result<(), std::io::Error> {
    let fq_iter = FqReader::new(path)?;
    let mut output = Output::new();
    if exclude_masked {
        for record in fq_iter.records() {
            match record {
                Ok(read) => {
                    let result = cal_unmasked_seq(&read);
                    print(&mut output, read.id(), read.seq().len(), &result)?;
                }
                Err(e) => eprintln!("Error read fASTQ: {}", e),
            }
        }
    } else {
        for record in fq_iter.records() {
            match record {
                Ok(read) => {
                    let result = cal_all_seq(&read);
                    print(&mut output, read.id(), read.seq().len(), &result)?;
                }
                Err(e) => eprintln!("Error read fASTQ: {}", e),
            }
        }
    }
    Ok(())
}
/// Parses FASTA file and compute the statistic w/o masked sequences.
/// Outputs the results to [`std::io::stdout()`].
///The output columns are:
/// - `#seq`: Number of reads
/// - `#bases`: Number of bases
/// - `#A`, `#C`, `#G`, `#T`: Number of each nucleotide
/// - `#2`: Number of `R`, `Y`, `S`, `W`, `K`, `M`
/// - `#3`: Number of `B`, `D`, `H`, `V`
/// - `#4`: Number of `N`
/// - `CG`: Number of `CG` on the template strand
/// - `GC`: Number of `GC` on the template strand
///
/// # Arguments
///
/// * `path` - FASTA path
/// * `exclude_masked` - If true, masked sequences (e.g., lowercases or char other then IUPAC code) will be excluded from the output.
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn calc_fa_comp_wo_bed(path: &str, exclude_masked: bool) -> Result<(), std::io::Error> {
    let fa_iter = FaReader::new(path)?;
    let mut output = Output::new();
    if exclude_masked {
        for record in fa_iter.records() {
            match record {
                Ok(read) => {
                    let result = cal_unmasked_seq(&read);
                    print(&mut output, read.id(), read.seq().len(), &result)?;
                }
                Err(e) => eprintln!("Error read fASTA: {}", e),
            }
        }
    } else {
        for record in fa_iter.records() {
            match record {
                Ok(read) => {
                    let result = cal_all_seq(&read);
                    print(&mut output, read.id(), read.seq().len(), &result)?;
                }
                Err(e) => eprintln!("Error read fASTA: {}", e),
            }
        }
    }
    Ok(())
}
/// Parses FASTA file and compute the statistic w/o masked sequences with a BED file.
/// Outputs the results to [`std::io::stdout()`].
///The output columns are:
/// - `#seq`: Number of reads
/// - `#bases`: Number of bases
/// - `#A`, `#C`, `#G`, `#T`: Number of each nucleotide
/// - `#2`: Number of `R`, `Y`, `S`, `W`, `K`, `M`
/// - `#3`: Number of `B`, `D`, `H`, `V`
/// - `#4`: Number of `N`
/// - `CG`: Number of `CG` on the template strand
/// - `GC`: Number of `GC` on the template strand
///
/// # Arguments
///
/// * `path` - FASTA path
/// * `bed` - BED path
/// * `exclude_masked` - If true, masked sequences (e.g., lowercases or char other then IUPAC code) will be excluded from the output.
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn calc_fq_comp_with_bed(
    path: &str,
    bed: &str,
    exclude_masked: bool,
) -> Result<(), std::io::Error> {
    let fq_iter = FqReader::new(path)?;
    let mut output = Output::new();
    if exclude_masked {
        let bedmap = BedMap::from(bed)?;
        for record in fq_iter.records() {
            match record {
                Ok(read) => {
                    if let Some(result) = cal_unmasked_seq_with_bed(&read, &bedmap) {
                        print(&mut output, read.id(), read.seq().len(), &result)?;
                    }
                }
                Err(e) => eprintln!("Error read fASTQ: {}", e),
            }
        }
    } else {
        let bedmap = BedMap::from(bed)?;
        for record in fq_iter.records() {
            match record {
                Ok(read) => {
                    if let Some(result) = cal_all_seq_with_bed(&read, &bedmap) {
                        print(&mut output, read.id(), read.seq().len(), &result)?;
                    }
                }
                Err(e) => eprintln!("Error read fASTQ: {}", e),
            }
        }
    }
    Ok(())
}
/// Parses FASTQ file and compute the statistic w/o masked sequences with a BED file.
/// Outputs the results to [`std::io::stdout()`].
///The output columns are:
/// - `#seq`: Number of reads
/// - `#bases`: Number of bases
/// - `#A`, `#C`, `#G`, `#T`: Number of each nucleotide
/// - `#2`: Number of `R`, `Y`, `S`, `W`, `K`, `M`
/// - `#3`: Number of `B`, `D`, `H`, `V`
/// - `#4`: Number of `N`
/// - `CG`: Number of `CG` on the template strand
/// - `GC`: Number of `GC` on the template strand
///
/// # Arguments
///
/// * `path` - FASTQ path
/// * `bed` - BED path
/// * `exclude_masked` - If true, masked sequences (e.g., lowercases or char other then IUPAC code) will be excluded from the output.
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn calc_fa_comp_with_bed(
    path: &str,
    bed: &str,
    exclude_masked: bool,
) -> Result<(), std::io::Error> {
    let fa_iter = FaReader::new(path)?;
    let mut output = Output::new();
    if exclude_masked {
        let bedmap = BedMap::from(bed)?;
        for record in fa_iter.records() {
            match record {
                Ok(read) => {
                    if let Some(result) = cal_unmasked_seq_with_bed(&read, &bedmap) {
                        print(&mut output, read.id(), read.seq().len(), &result)?;
                    }
                }
                Err(e) => eprintln!("Error read fASTA: {}", e),
            }
        }
    } else {
        let bedmap = BedMap::from(bed)?;
        for record in fa_iter.records() {
            match record {
                Ok(read) => {
                    if let Some(result) = cal_all_seq_with_bed(&read, &bedmap) {
                        print(&mut output, read.id(), read.seq().len(), &result)?;
                    }
                }
                Err(e) => eprintln!("Error read fASTA: {}", e),
            }
        }
    }
    Ok(())
}
fn cal_all_seq<T: RecordType>(read: &T) -> [usize; 9] {
    let mut count: [usize; 23] = [0; 23];
    SeqComp::count_all_nc(&mut count, read.seq(), 0, read.seq().len());
    SeqComp::get_all_result(&count)
}
fn cal_unmasked_seq<T: RecordType>(read: &T) -> [usize; 9] {
    let mut count: [usize; 23] = [0; 23];
    SeqComp::count_unmasked_nc(&mut count, read.seq(), 0, read.seq().len());
    SeqComp::get_unmasked_result(&count)
}
fn cal_all_seq_with_bed<T: RecordType>(read: &T, bedmap: &BedMap) -> Option<[usize; 9]> {
    let mut count: [usize; 23] = [0; 23];
    if let Some(bedvec) = bedmap.get(read.id()) {
        bedvec.iter().for_each(|pos| {
            SeqComp::count_all_nc(&mut count, read.seq(), pos.0, pos.1);
        });
        let result = SeqComp::get_all_result(&count);
        return Some(result);
    }
    None
}
fn cal_unmasked_seq_with_bed<T: RecordType>(read: &T, bedmap: &BedMap) -> Option<[usize; 9]> {
    let mut count: [usize; 23] = [0; 23];
    if let Some(bedvec) = bedmap.get(read.id()) {
        bedvec.iter().for_each(|pos| {
            SeqComp::count_unmasked_nc(&mut count, read.seq(), pos.0, pos.1);
        });
        let result = SeqComp::get_unmasked_result(&count);
        return Some(result);
    }
    None
}
fn print(
    output: &mut Output,
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
