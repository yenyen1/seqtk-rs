use crate::dna::DNA;

pub struct SequenceRead {
    name: String,
    seq: Vec<DNA>,
    qual: String,
    ascii_bases: u8,
}

impl SequenceRead {
    pub fn empty_read() -> SequenceRead {
        SequenceRead {
            name: String::new(),
            seq: Vec::new(),
            qual: String::new(),
            ascii_bases: 0,
        }
    }
    pub fn set_name(&mut self, name: &str) {
        self.name = name.to_string();
    }
    pub fn set_seq(&mut self, seq_string: &str) {
        self.seq = seq_string
            .chars()
            .map(|c| match c {
                'A' => DNA::A,
                'a' => DNA::A,
                'T' => DNA::T,
                't' => DNA::T,
                'C' => DNA::C,
                'c' => DNA::C,
                'G' => DNA::G,
                'g' => DNA::G,
                _ => DNA::N,
            })
            .collect();
    }
    pub fn set_qual(&mut self, qual_string: &str, ascii_base: u8) {
        self.ascii_bases = ascii_base;
        self.qual = qual_string.to_string();
        // self.qual = qual_string.chars().map(|c| c as u8 - ascii_base).collect();
    }
    pub fn get_seq(&self) -> &Vec<DNA> {
        &self.seq
    }
    pub fn get_name(&self) -> &str {
        self.name.split_whitespace().next().unwrap_or("")
    }
    pub fn get_q_qual_score_u8(&self) -> Vec<u8> {
        self.qual
            .chars()
            .map(|c| c as u8 - self.ascii_bases)
            .collect()
    }

    pub fn get_p_error_score_f32(&self) -> Vec<f32> {
        self.qual
            .chars()
            .map(|c| 10.0f32.powf(-0.1 * (c as u8 - self.ascii_bases) as f32))
            .collect()
    }

    pub fn get_qaul_char(&self) -> Vec<char> {
        self.qual.chars().collect()
    }

    pub fn get_seq_len(&self) -> usize {
        self.seq.len()
    }
}
