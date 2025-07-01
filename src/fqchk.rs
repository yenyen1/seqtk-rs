use crate::io_utils::{FqReader, Output};
use crate::stats;
use rayon::prelude::*;
use std::fmt::Write;

/// Parses FASTQ data without quality threshold and computes per-position statistics.
/// Outputs the results to [`std::io::stdout()`].
///
/// The output columns are:
/// - `POS`: Position in the read
/// - `#bases`: Number of bases at this position
/// - `%A`, `%C`, `%G`, `%T`, `%N`: Percentage of each nucleotide
/// - `avgQ`: Average quality score `(Q₁ + Q₂ + ... + Qₙ) / N`
/// - `errQ`: Estimated error rate `-10 * log₁₀((P₁ + P₂ + ... + Pₙ) / N)`
/// - `%Qx`: Percentage of each quality score
///
/// # Arguments
///
/// * `path` - FASTQ path
/// * `asciibase` - Quality scores equal to the score plus a base offset asciibase.
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn get_fqchk_result_wo_qual_threshold(
    path: &str,
    asciibase: usize,
) -> Result<(), std::io::Error> {
    let (maxlen, qualset) = get_maxlen_and_qualset(path)?;
    cal_seq_all(path, maxlen, &qualset, asciibase)?;
    Ok(())
}

/// Parses FASTQ data with quality threshold and computes per-position statistics.
/// Outputs the results to [`std::io::stdout()`].
///
/// The output columns are:
/// - `POS`: Position in the read
/// - `#bases`: Number of bases at this position
/// - `%A`, `%C`, `%G`, `%T`, `%N`: Percentage of each nucleotide
/// - `avgQ`: Average quality score `(Q₁ + Q₂ + ... + Qₙ) / N`
/// - `errQ`: Estimated error rate `-10 * log₁₀((P₁ + P₂ + ... + Pₙ) / N)`
/// - `%low`, `%high`: Percentage of the nucleotide that the quality scores below or above the threshold, respectively.
///
/// # Arguments
///
/// * `path` - FASTQ path
/// * `q_plus_ascii` - The sum of quality threshold and asciibase.
/// * `asciibase` - Quality score equal to the score plus a base offset asciibase.
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn get_fqchk_with_qual_threshold(
    path: &str,
    q_plus_ascii: u8,
    asciibase: usize,
) -> Result<(), std::io::Error> {
    let (maxlen, qualset) = get_maxlen_and_qualset(path)?;
    cal_seq_with_q(path, maxlen, &qualset, asciibase, q_plus_ascii)?;
    Ok(())
}

fn cal_seq_all(
    path: &str,
    maxlen: usize,
    qual_set: &[usize],
    asciibases: usize,
) -> Result<(), std::io::Error> {
    let mut seq_count_mat: Vec<[usize; 256]> = vec![[0; 256]; maxlen];
    let mut seq_pos_sum: Vec<usize> = vec![0; maxlen];
    let mut seq_all: [usize; 256] = [0; 256];
    let mut qual_count_mat: Vec<[usize; 256]> = vec![[0; 256]; maxlen];
    let mut qual_all: [usize; 256] = [0; 256];

    let fq = FqReader::new(path)?;
    for record in fq.records() {
        match record {
            Ok(read) => {
                for (i, (&b, &q)) in read.seq().iter().zip(read.qual()).enumerate() {
                    let ub = b as usize;
                    seq_count_mat[i][ub] += 1;
                    seq_all[ub] += 1;
                    seq_pos_sum[i] += 1;
                    let uq = q as usize;
                    qual_count_mat[i][uq] += 1;
                    qual_all[uq] += 1;
                }
            }
            Err(e) => eprintln!("Error read FASTQ: {}", e),
        }
    }

    let out = std::io::stdout();
    let mut output = Output::new(out);

    let column = format!(
        "{}\t{}\n",
        "POS\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ",
        get_qual_cols(qual_set, asciibases)
    );
    output.write(column)?;

    let mut total = seq_all.iter().sum();
    let mut total_f64 = total as f64;
    let mut result = format!(
        "all\t{}\t{}\t{}\n",
        get_seq_result(total, &seq_all),
        get_avg_err(total_f64, &qual_all, qual_set, asciibases),
        get_qual_result(total_f64, &qual_all, qual_set),
    );
    output.write(result)?;
    for i in 0..maxlen {
        total = seq_pos_sum[i];
        total_f64 = total as f64;
        result = format!(
            "{}\t{}\t{}\t{}\n",
            i + 1,
            get_seq_result(total, &seq_count_mat[i]),
            get_avg_err(total_f64, &qual_count_mat[i], qual_set, asciibases),
            get_qual_result(total_f64, &qual_count_mat[i], qual_set),
        );
        output.write(result)?;
    }
    Ok(())
}
fn cal_seq_with_q(
    path: &str,
    maxlen: usize,
    qual_set: &[usize],
    asciibases: usize,
    q_plus_ascii: u8,
) -> Result<(), std::io::Error> {
    let mut seq_count_mat: Vec<[usize; 256]> = vec![[0; 256]; maxlen];
    let mut seq_all: [usize; 256] = [0; 256];
    let mut qual_count_mat: Vec<[usize; 256]> = vec![[0; 256]; maxlen];
    let mut qual_all: [usize; 256] = [0; 256];
    let mut qual_q_count: Vec<[usize; 2]> = vec![[0; 2]; maxlen]; // low, high
    let mut qual_q_count_all: [usize; 2] = [0; 2];

    let fq = FqReader::new(path)?;
    for record in fq.records() {
        match record {
            Ok(read) => {
                for (i, (&b, &q)) in read.seq().iter().zip(read.qual()).enumerate() {
                    let ub = b as usize;
                    seq_count_mat[i][ub] += 1;
                    seq_all[ub] += 1;
                    let uq = q as usize;
                    qual_count_mat[i][uq] += 1;
                    qual_all[uq] += 1;
                    if q >= q_plus_ascii {
                        qual_q_count[i][1] += 1;
                        qual_q_count_all[1] += 1;
                    } else {
                        qual_q_count[i][0] += 1;
                        qual_q_count_all[0] += 1;
                    }
                }
            }
            Err(e) => eprintln!("Error read fASTQ: {}", e),
        }
    }

    let out = std::io::stdout();
    let mut output = Output::new(out);
    let column = format!(
        "{}\t{}\n",
        "POS\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ", "%low\t%high",
    );
    output.write(column)?;
    let mut total = qual_q_count_all[0] + qual_q_count_all[1];
    let mut total_f64 = total as f64;
    let mut result = format!(
        "all\t{}\t{}\t{}\n",
        get_seq_result(total, &seq_all),
        get_avg_err(total_f64, &qual_all, qual_set, asciibases),
        get_qual_with_q_result(total_f64, &qual_q_count_all),
    );
    output.write(result)?;
    for i in 0..maxlen {
        total = qual_q_count[i][0] + qual_q_count[i][1];
        total_f64 = total as f64;
        result = format!(
            "{}\t{}\t{}\t{}\n",
            i + 1,
            get_seq_result(total, &seq_count_mat[i]),
            get_avg_err(total_f64, &qual_count_mat[i], qual_set, asciibases),
            get_qual_with_q_result(total_f64, &qual_q_count[i]),
        );
        output.write(result)?;
    }
    Ok(())
}
fn get_maxlen_and_qualset(path: &str) -> Result<(usize, Vec<usize>), std::io::Error> {
    let fq = FqReader::new(path)?;
    let mut maxlen: usize = 0;
    let mut qual_set: [bool; 256] = [false; 256];
    for record in fq.records() {
        match record {
            Ok(read) => {
                let len = read.seq().len();
                if len > maxlen {
                    maxlen = len;
                }

                for &qual in read.qual() {
                    qual_set[qual as usize] = true;
                }
            }
            Err(e) => eprintln!("Error read FASTQ: {}", e),
        }
    }
    let uniq_qset: Vec<usize> = qual_set
        .iter()
        .enumerate()
        .filter_map(|(i, &b)| if b { Some(i) } else { None })
        .collect();
    Ok((maxlen, uniq_qset))
}

/// [Note] Some tools treat Q < 3 as Q = 3. I don't do that.
/// Q = 0 leads to P = 1.0 ; Q = 1 → P = 0.794 ; Q = 2 → P = 0.630.
/// These small Qs significantly affect and skew the result of errQ.
/// Therefore, they treat Q < 3 as Q = 3.
fn get_avg_err(
    total: f64,
    qual_count: &[usize; 256],
    qual_set: &[usize],
    asciibases: usize,
) -> String {
    let sum: f64 = qual_set
        .par_iter()
        .map(|&q| ((q - asciibases) as f64) * (qual_count[q] as f64))
        .sum();
    let avg_q = sum / total;
    let sum: f64 = qual_set
        .par_iter()
        .map(|&q| stats::convert_q_score_to_p_err((q - asciibases) as f64) * (qual_count[q] as f64))
        .sum();
    let err_q = stats::convert_p_err_to_q_score(sum / total);

    format!("{:.1}\t{:.1}", avg_q, f64::abs(err_q))
}
/// Output: %Qx
fn get_qual_result(total: f64, qual_count: &[usize; 256], qual_set: &[usize]) -> String {
    let mut result = String::new();
    for (i, &q) in qual_set.iter().enumerate() {
        if i > 0 {
            result.push('\t');
        }
        write!(&mut result, "{:.1}", qual_count[q] as f64 * 100.0 / total).unwrap();
    }
    result
}
/// Output: %low, %high
fn get_qual_with_q_result(total_f64: f64, qual_count: &[usize; 2]) -> String {
    let result = format!(
        "{:.1}\t{:.1}",
        qual_count[0] as f64 * 100.0 / total_f64,
        qual_count[1] as f64 * 100.0 / total_f64
    );
    result
}
fn get_seq_result(total: usize, seq_count: &[usize; 256]) -> String {
    let total_f64 = total as f64;
    let result = format!(
        "{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}",
        total,
        100.0 * (seq_count[b'A' as usize] + seq_count[b'a' as usize]) as f64 / total_f64,
        100.0 * (seq_count[b'C' as usize] + seq_count[b'c' as usize]) as f64 / total_f64,
        100.0 * (seq_count[b'G' as usize] + seq_count[b'g' as usize]) as f64 / total_f64,
        100.0 * (seq_count[b'T' as usize] + seq_count[b't' as usize]) as f64 / total_f64,
        100.0 * (seq_count[b'N' as usize] + seq_count[b'n' as usize]) as f64 / total_f64,
    );
    result
}
fn get_qual_cols(qual_set: &[usize], asciibases: usize) -> String {
    let mut result = String::new();
    for (i, &q) in qual_set.iter().enumerate() {
        if i > 0 {
            result.push('\t');
        }
        write!(&mut result, "%Q{}", q - asciibases).unwrap();
    }
    result
}
