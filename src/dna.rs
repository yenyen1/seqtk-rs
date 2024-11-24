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
    // pub fn all_variants() -> Vec<DNA> {
    //     vec![DNA::A, DNA::C, DNA::G, DNA::T, DNA::N]
    // }
    pub fn all_chars() -> Vec<char> {
        vec!['A', 'C', 'G', 'T', 'N']
    }
    // pub fn get_dna_idx(dna: &DNA) -> usize {
    //     match *dna {
    //         DNA::A => 0,
    //         DNA::C => 1,
    //         DNA::G => 2,
    //         DNA::T => 3,
    //         DNA::N => 4,
    //     }
    // }
    pub fn get_char_idx(c: &char) -> usize {
        match *c {
            'A' | 'a' => 0,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            _ => 4,
        }
    }

    // pub fn convert_char_to_dna(seq_char: &char) -> DNA {
    //     match seq_char {
    //         'A' => DNA::A,
    //         'a' => DNA::A,
    //         'T' => DNA::T,
    //         't' => DNA::T,
    //         'C' => DNA::C,
    //         'c' => DNA::C,
    //         'G' => DNA::G,
    //         'g' => DNA::G,
    //         _ => DNA::N,
    //     }
    // }
    // pub fn convert_dna_to_char(dna: &DNA) -> char {
    //     match *dna {
    //         DNA::A => 'A',
    //         DNA::C => 'C',
    //         DNA::G => 'G',
    //         DNA::T => 'T',
    //         DNA::N => 'N',
    //     }
    // }
}
