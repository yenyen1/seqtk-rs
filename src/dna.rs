use std::collections::HashMap;

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
    pub fn all_chars() -> Vec<char> {
        DNA::all_variants()
            .iter()
            .map(|d| DNA::convert_dna_to_char(d))
            .collect()
    }
    pub fn get_dna(idx: usize) -> DNA {
        match idx {
            0 => DNA::A,
            1 => DNA::C,
            2 => DNA::G,
            3 => DNA::T,
            _ => DNA::N,
        }
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
    pub fn convert_dna_to_char(dna: &DNA) -> char {
        match *dna {
            DNA::A => 'A',
            DNA::C => 'C',
            DNA::G => 'G',
            DNA::T => 'T',
            DNA::N => 'N',
        }
    }
}
