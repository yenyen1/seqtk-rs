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
    // A helper function to return a list of all variants
    pub fn all_variants() -> Vec<DNA> {
        vec![DNA::A, DNA::C, DNA::G, DNA::T, DNA::N]
    }

    pub fn convert_char_to_dna(seq_char: &char) -> DNA {
        match seq_char {
            'A' => DNA::A,
            'a' => DNA::A,
            'T' => DNA::T,
            't' => DNA::T,
            'C' => DNA::C,
            'c' => DNA::C,
            'G' => DNA::G,
            'g' => DNA::G,
            _ => DNA::N,
        }
    }
}
