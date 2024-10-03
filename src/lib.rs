pub fn abc(fq_path: &str, quality_value: Option<u32>) {
    let q = quality_value.unwrap_or(0);
    println!("abc: {} {}", fq_path, q);
}
