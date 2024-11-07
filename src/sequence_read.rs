use crate::dna::DNA;

pub struct SequenceRead {
    name: String,
    seq: String,
    qual: String,
    ascii_bases: u8,
}

impl SequenceRead {
    /// [Construct SequenceRead] 
    /// Create empty SequenceRead 
    pub fn empty_read(ascii_base: u8) -> SequenceRead {
        SequenceRead {
            name: String::new(),
            seq: String::new(),
            qual: String::new(),
            ascii_bases: ascii_base,
        }
    }
    /// Update SequenceRead with input String
    pub fn set_read_name(&mut self, name: &str) {
        self.name = name.to_string();
    }
    pub fn set_read_seq(&mut self, seq_string: &str) {
        self.seq = seq_string.to_string();
    }
    pub fn set_read_qual_str(&mut self, qual_string: &str) {
        self.qual = qual_string.to_string();
    }
    /// Get SequenceRead data
    pub fn get_qual_chars(&self) -> Vec<char> {
        self.qual.chars().collect()
    }
    pub fn get_seq_length(&self) -> usize {
        self.seq.len()
    }
    pub fn get_seq(&self) -> Vec<DNA> {
        self.seq
                .chars()
                .map(|c| DNA::convert_char_to_dna(&c))
                .collect()
    }
    
    
    // pub fn get_name(&self) -> &str {
    //     self.name.split_whitespace().next().unwrap_or("")
    // }
    pub fn get_q_qual_score_u8(&self) -> Vec<u8> {
        self.qual
            .chars()
            .map(|c| c as u8 - self.ascii_bases)
            .collect()
    }

    // pub fn get_p_error_score_f64(&self) -> Vec<f64> {
    //     self.qual
    //         .chars()
    //         .map(|c| 10.0f64.powf(-0.1 * (c as u8 - self.ascii_bases) as f64))
    //         .collect()
    // }
}
