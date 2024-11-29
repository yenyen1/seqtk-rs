pub fn ordered_dna_str() -> String {
    "\tA\tC\tG\tT\tN".to_string()
}
pub fn get_char_idx(c: char) -> usize {
    match c {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        _ => 4,
    }
}
