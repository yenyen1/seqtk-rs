use crate::dna::DNA;

pub struct SequenceRead {
    name: String,
    seq: String,
    qual: String,
    ascii_bases: u8,
}

impl SequenceRead {
    pub fn empty_read() -> SequenceRead {
        SequenceRead {
            name: String::new(),
            seq: String::new(),
            qual: String::new(),
            ascii_bases: 0,
        }
    }
    pub fn set_name(&mut self, name: &str) {
        self.name = name.to_string();
    }
    pub fn set_seq(&mut self, seq_string: &str) {
        self.seq = seq_string.to_string();
    }
    pub fn set_qual(&mut self, qual_string: &str, ascii_base: u8) {
        self.ascii_bases = ascii_base;
        self.qual = qual_string.to_string();
    }
    pub fn get_seq_with_qual_value_threshold(&self, qual_threshold: u8) -> Vec<DNA> {
        match qual_threshold {
            0 => self
                .seq
                .chars()
                .map(|c| DNA::convert_char_to_dna(&c))
                .collect(),
            _ => {
                self.seq
                    .chars()
                    .zip(self.qual.chars())
                    .map(|(c, q)| {
                        self.convert_char_to_dna_with_qual_threshold(&c, &q, &qual_threshold)
                    })
                    .collect()
            }
        }
    }
    fn convert_char_to_dna_with_qual_threshold(
        &self,
        &seq_char: &char,
        &qual_char: &char,
        &qual_threshold: &u8,
    ) -> DNA {
        if (qual_char as u8 - self.ascii_bases) >= qual_threshold {
            DNA::convert_char_to_dna(&seq_char)
        } else {
            DNA::N
        }
    }

    // pub fn get_name(&self) -> &str {
    //     self.name.split_whitespace().next().unwrap_or("")
    // }
    // pub fn get_q_qual_score_u8(&self) -> Vec<u8> {
    //     self.qual
    //         .chars()
    //         .map(|c| c as u8 - self.ascii_bases)
    //         .collect()
    // }

    // pub fn get_p_error_score_f64(&self) -> Vec<f64> {
    //     self.qual
    //         .chars()
    //         .map(|c| 10.0f64.powf(-0.1 * (c as u8 - self.ascii_bases) as f64))
    //         .collect()
    // }

    pub fn get_qaul_char(&self) -> Vec<char> {
        self.qual.chars().collect()
    }

    pub fn get_seq_len(&self) -> usize {
        self.seq.len()
    }
}
