#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Hash, Eq, PartialEq)]

pub enum DNA {
    A,
    C,
    G,
    T,
    N,
}

impl DNA {
    pub fn ordered_dna_str() -> String {
        format!("\tA\tC\tG\tT\tN")
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
}
