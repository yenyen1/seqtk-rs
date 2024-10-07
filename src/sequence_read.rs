use crate::dna::DNA;

pub struct SequenceRead {
    name: String,
    seq: Vec<DNA>,
    qual: Vec<u8>,
}

impl SequenceRead {
    pub fn empty_read() -> SequenceRead {
        SequenceRead {
            name: String::new(),
            seq: Vec::new(),
            qual: Vec::new(),
        }
    }
    pub fn set_name(&mut self, name: &str) {
        self.name = name.split_whitespace().next().unwrap_or("").to_string();
    }
    pub fn set_seq(&mut self, seq_string: &str) {
        self.seq = seq_string
            .chars()
            .map(|c| match c {
                'A' => DNA::A(true),
                'a' => DNA::A(false),
                'T' => DNA::T(true),
                't' => DNA::T(false),
                'C' => DNA::C(true),
                'c' => DNA::C(false),
                'G' => DNA::G(true),
                'g' => DNA::G(false),
                _ => DNA::N,
            })
            .collect();
    }
    pub fn set_qual(&mut self, qual_string: &str, ascii_base: u8) {
        match ascii_base {
            33 | 64 => {
                self.qual = qual_string.chars().map(|c| c as u8 - ascii_base).collect();
            }
            _ => {
                println!("ascii base not support");
            }
        }
    }
    pub fn get_seq(&self) -> &Vec<DNA> {
        &self.seq
    }
    pub fn get_name(&self) -> &str {
        self.name.as_str()
    }
    pub fn get_qual(&self) -> &Vec<u8> {
        &self.qual
    }
    pub fn get_seq_len(&self) -> usize {
        self.seq.len()
    }
}
