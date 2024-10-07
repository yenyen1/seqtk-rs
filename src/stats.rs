pub fn median(sorted_arr: &[usize]) -> f64 {
    let len_a = sorted_arr.len();
    match len_a {
        0 => f64::NAN,
        _ => {
            let mid = len_a / 2;
            match len_a % 2 {
                1 => sorted_arr[mid] as f64,
                _ => ((sorted_arr[mid - 1] as f64) + (sorted_arr[mid] as f64)) / 2.0,
            }
        }
    }
}

pub fn average(arr: &[usize]) -> f64 {
    let len_a = arr.len();
    match len_a {
        0 => f64::NAN,
        _ => arr.iter().sum::<usize>() as f64 / len_a as f64,
    }
}
pub fn n50(sorted_arr: &[usize]) -> usize {
    let len_a = sorted_arr.len();
    match len_a {
        0 => 0,
        _ => {
            let half_sum = sorted_arr.iter().sum::<usize>() / 2;
            let mut cur_acc_bases = 0;
            for cur_read_size in sorted_arr.iter().rev() {
                cur_acc_bases += cur_read_size;
                if cur_acc_bases >= half_sum {
                    return *cur_read_size;
                }
            }
            panic!("Should not hit this line");
        }
    }
}
