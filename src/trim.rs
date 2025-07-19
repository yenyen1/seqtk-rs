use crate::stats::Q2PConverter;
pub fn trimfq() {}

/// return not_qualified, length
fn trim_read_by_error_rate(
    asciibase: u8,
    qual: &[u8],
    ethershold: f64,
    minlen: usize,
) -> (bool, usize) {
    let mut err_sum = vec![0.0; qual.len()];
    let qp_converter = Q2PConverter::new(asciibase);
    // let mut p = stats::convert_q_score_to_p_err((qual[0]-asciibase) as f64);
    err_sum[0] = qp_converter.get_prob(qual[0]);
    (1..qual.len()).for_each(|i| {
        // p = stats::convert_q_score_to_p_err((qual[i]-asciibase) as f64);
        err_sum[i] = err_sum[i - 1] + qp_converter.get_prob(qual[i]);
    });

    let avg_sum: Vec<f64> = err_sum
        .iter()
        .enumerate()
        .map(|(i, p)| p / (i + 1) as f64)
        .collect();
    println!("thres: {} sum: {:?}", ethershold, avg_sum);

    for i in (0..err_sum.len()).rev() {
        if err_sum[i] <= ethershold * (i + 1) as f64 {
            return (false, i + 1);
        }
        if i < minlen {
            return (true, minlen);
        }
    }
    (true, minlen)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats;

    #[test]
    fn test_trim_read() {
        // let qual: [u8; 20] = [40, 40, 40, 0, 0, 40, 40, 40, 40, 40, 0, 0, 40, 40, 40, 40, 40, 40, 40, 40];
        // let (_, size) = trim_read_by_error_rate(0, &qual, 0.05*20.0, 0);
        // assert_eq!(size, 3);

        let qual = b"...&...++.&.++...&...''..'.&..+..++++.++";
        // let qual_t: Vec<f64> = qual.iter().map(|q| (q-33) as f64).collect();
        // let p_err: Vec<f64> = qual_t.iter().enumerate().map(|(i, &q)| stats::convert_q_score_to_p_err(q)/(i+1)as f64 ).collect();
        // println!("{:?}", p_err);
        let (_, size) = trim_read_by_error_rate(33, qual, 0.1, 0);
        // assert_eq!(size, 38);
    }
}
