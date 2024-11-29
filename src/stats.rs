use std::iter::Sum;
use std::ops::Add;
// use std::iter::Iterator::Sum;

pub fn median_sorted_arr<T>(sorted_arr: &[T]) -> f64
where
    T: Copy + Into<f64>,
{
    let len_a = sorted_arr.len();
    match len_a {
        0 => f64::NAN,
        _ => {
            let mid = len_a / 2;
            match len_a % 2 {
                1 => sorted_arr[mid].into(),
                _ => (sorted_arr[mid - 1].into() + sorted_arr[mid].into()) / 2.0,
            }
        }
    }
}

pub fn average<T>(arr: &[T]) -> f64
where
    T: Copy + Add<Output = T> + Sum + Into<f64>,
{
    let len_a = arr.len();
    if len_a == 0 {
        return f64::NAN;
    }
    let sum: f64 = arr.iter().copied().map(Into::into).sum();

    sum / (len_a as f64)
}

pub fn n50_sorted_arr<T>(sorted_arr: &[T]) -> Option<T>
where
    T: Copy + Add<Output = T> + Sum + Into<u64>,
{
    let len_a = sorted_arr.len();
    if len_a == 0 {
        return None;
    }

    let total_sum: u64 = sorted_arr.iter().copied().map(Into::into).sum();
    let half_sum = total_sum / 2;
    let mut cur_acc_bases: u64 = 0;
    for &cur_read_size in sorted_arr.iter().rev() {
        cur_acc_bases += cur_read_size.into();
        if cur_acc_bases >= half_sum {
            return Some(cur_read_size);
        }
    }

    // If the function reaches here, it indicates an unexpected state
    None
}
pub fn convert_p_err_to_q_score(p_err: f64) -> f64 {
    -10.0 * p_err.log10()
}
pub fn convert_q_score_to_p_err(q_score: f64) -> f64 {
    10.0f64.powf(-0.1 * q_score)
}
