use strum_macros::EnumIter;

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, EnumIter, Hash, Clone, Eq, PartialEq)]
pub enum DNA {
    A,
    T,
    C,
    G,
    N,
}

impl DNA {
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
