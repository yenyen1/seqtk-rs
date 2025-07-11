use crate::io_utils::{FaReader, FqReader, Output};
use rayon::slice::ParallelSliceMut;

/// Parses FASTQ file and computes read statistics.
/// Outputs the results to [`std::io::stdout()`].
///The output columns are:
/// - `#seq`: Number of reads
/// - `#bases`: Total number of bases
/// - `avg_size`: Average read length
/// - `mini_size`: Minimum read length
/// - `med_size`: Median read length
/// - `max_size`: Maximum read length
/// - `N50`: N50 read length
///
/// # Arguments
///
/// * `path` - FASTQ path
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn calc_fq_size(path: &str) -> Result<(), std::io::Error> {
    let fq_iter = FqReader::new(path)?;
    let mut seq_len: Vec<usize> = Vec::new();
    for record in fq_iter.records() {
        match record {
            Ok(read) => {
                seq_len.push(read.seq().len());
            }
            Err(e) => eprintln!("Error read fASTQ: {}", e),
        }
    }
    seq_len.par_sort_unstable();
    let result = get_result_str(&seq_len);
    let mut output = Output::new();

    output.write(result)?;
    Ok(())
}
/// Parses FASTA file and computes sequence statistics.
/// Outputs the results to [`std::io::stdout()`].
///The output columns are:
/// - `#seq`: Number of sequences
/// - `#bases`: Total number of bases
/// - `avg_size`: Average sequence length
/// - `mini_size`: Minimum sequence length
/// - `med_size`: Median sequence length
/// - `max_size`: Maximum sequence length
/// - `N50`: N50 sequence length
///
/// # Arguments
///
/// * `path` - FASTA path
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn calc_fa_size(path: &str) -> Result<(), std::io::Error> {
    let fa_iter = FaReader::new(path)?;
    let mut seq_len: Vec<usize> = Vec::new();
    for record in fa_iter.records() {
        match record {
            Ok(read) => {
                seq_len.push(read.seq().len());
            }
            Err(e) => eprintln!("Error read fASTA: {}", e),
        }
    }
    seq_len.par_sort_unstable();
    let result = get_result_str(&seq_len);
    let mut output = Output::new();
    output.write(result)?;
    Ok(())
}

fn get_result_str(sorted_seq_len: &[usize]) -> String {
    // #seq, #bases, avg_size, min_size, med_size, max_size, N50
    let sum: usize = sorted_seq_len.iter().sum();
    let size = sorted_seq_len.len();
    let median = match size {
        0 => f64::NAN,
        _ => {
            let mid = size / 2;
            match size % 2 {
                1 => sorted_seq_len[mid] as f64,
                _ => (sorted_seq_len[mid - 1] + sorted_seq_len[mid]) as f64 / 2.0,
            }
        }
    };
    let mut n50: usize = 0;
    let half: usize = sum / 2;
    let mut acc: usize = 0;
    for &cur in sorted_seq_len.iter().rev() {
        acc += cur;
        if acc >= half {
            n50 = cur;
            break;
        }
    }
    let result = format!(
        "{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
        sorted_seq_len.len(),
        sum,
        sum as f64 / size as f64,
        sorted_seq_len.first().unwrap(),
        median,
        sorted_seq_len.last().unwrap(),
        n50
    );
    result
}
