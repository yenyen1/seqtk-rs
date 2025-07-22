use crate::io_utils::{FqReader, FxWriter};
/// Trims low-quality bases from a FASTQ data based on a quality threshold `Q`.
/// Outputs the FASTQ to [`std::io::stdout()`].
///
/// The algorithm:
///
/// (1) Scan the read from 5' to 3' to find the first base with `quality >= Q`.
///     This is the starting position of the trimmed read.
/// (2) Compute a running score from 5' to 3'
///        `score(i) = score(i-1) + quality(i) - Q`
///      where the initial score is `score(start) = quality(start) - Q`.
/// (3) Find out the maximun score. This is the end of the trimmed read.
///
/// Note:
/// (1) If all base quality are less than `Q`, the read is discarded.
/// (2) If the trimmed read is shorter than `min_len`, we first extend it from the 3' end. If still too short, extend form the 5' end until `min_len` is reached or no more bases are available.
///
///
/// # Arguments
///
/// * `path` - FASTQ path
/// * `q_plus_ascii` - The sum of quality threshold and asciibase.
/// * `minilen` - The minimum length of read.
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
///
/// # Notes
/// Quality trimming is no longer necessary in most modern sequencing pipelines.
/// Its usefulness depends on the sequencing technology and the goals of your downstream analysis.
///
/// For example:
/// * Illumina reads may benefit from light trimming low-quality tails at 3' end.
/// * Long-read technologies are generally not necessary to do quality trimming, because it may remove informative regions and reduce read length unnecessarily.
///
pub fn trimfq(fq_path: &str, q_plus_ascii: u8, minlen: usize) -> Result<(), std::io::Error> {
    let reader = FqReader::new(fq_path)?;
    let mut writer = FxWriter::new(false);
    let (mut start, mut end, mut size): (usize, usize, usize);
    for record in reader.records() {
        let read = record.unwrap();
        (start, end) = trim_read_by_q(read.qual(), q_plus_ascii);

        // discard read if start == end (Q < q_threshold for all bases)
        if start < end {
            size = read.qual().len();

            if size < minlen {
                // write full length is size < minlen
                (start, end) = (0, size - 1);
            } else if minlen > (end - start + 1) {
                if size - start >= minlen {
                    // append 3' seq to reach the minlen
                    end = start + minlen - 1;
                } else {
                    // append 5' seq if the whole 3' end is not enough
                    end = size - 1;
                    start = size - minlen;
                }
            }
            writer.write(
                read.id(),
                &read.seq()[start..(end + 1)],
                read.desc(),
                &read.qual()[start..(end + 1)],
            )?;
        }
    }
    Ok(())
}

/// return start, end
fn trim_read_by_q(qual: &[u8], q_plus_ascii: u8) -> (usize, usize) {
    let mut start: usize = 0;
    let len = qual.len();
    for &q in qual {
        if q >= q_plus_ascii {
            break;
        } else {
            start += 1;
        }
    }
    if start == len {
        return (len, len);
    }

    let qthd_add_ascii_i = q_plus_ascii as isize;
    let mut err_sum: Vec<isize> = vec![0; len];

    err_sum[start] = qual[start] as isize - qthd_add_ascii_i;
    let (mut max_idx, mut max_val) = (start, err_sum[start]);

    ((start + 1)..len).for_each(|i| {
        err_sum[i] = err_sum[i - 1] + qual[i] as isize - qthd_add_ascii_i;
        if err_sum[i] >= max_val {
            max_idx = i;
            max_val = err_sum[i];
        }
    });

    // println!("[!] thres: {} sum: {:?}", qthd_add_ascii_i, &err_sum[start..]);

    (start, max_idx)
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_read() {
        let qual = b"&&&&&...&...++.&.++...&...''..'.&..+..++++.++";
        let (start, end) = trim_read_by_q(qual, 33 + 10);
        assert_eq!(end - start + 1, 40);
        let qual= b"?@<DDD;2A<><FIBHB?FCGAHHEBEHAFCBEFGGB@4CCA@?*?DH9?B?BDC/?81.8BFFEHG@>@G=DGCECEHEHF>>;?;?B=3@CC:(,(98',>5(82)<5>>223@4+4>@+5855>>>@:AA@>:43>@##########";

        let (start, end) = trim_read_by_q(qual, 33 + 20);
        assert_eq!(end - start + 1, 140);
    }
}
