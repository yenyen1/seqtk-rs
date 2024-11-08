use crate::dna::DNA;

pub struct SequenceRead {
    name: String,
    seq: String,
    qual: Vec<u8>,
    ascii_bases: u8,
}

impl SequenceRead {
    /// [Construct SequenceRead]
    /// Create empty SequenceRead
    pub fn empty_read(ascii_base: u8) -> SequenceRead {
        SequenceRead {
            name: String::new(),
            seq: String::new(),
            qual: Vec::<u8>::new(),
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
        self.qual = qual_string.chars().map(|c| c as u8 - self.ascii_bases).collect();
    }
    /// Get SequenceRead data
    pub fn get_seq_length(&self) -> usize {
        self.seq.len()
    }
    pub fn get_seq(&self) -> Vec<DNA> {
        self.seq
            .chars()
            .map(|c| DNA::convert_char_to_dna(&c))
            .collect()
    }
    pub fn get_seq_char(&self, idx: usize) -> DNA {
        let seq_char = self.seq.as_bytes()[idx] as char;
        DNA::convert_char_to_dna(&seq_char)
    }

    // pub fn get_name(&self) -> &str {
    //     self.name.split_whitespace().next().unwrap_or("")
    // }
    pub fn get_q_score_vec(&self) -> &Vec<u8> {
        &self.qual
    }
    pub fn get_p_err_vec(&self) -> Vec<f64> {
        self.qual.iter()
            .map(|c| SequenceRead::convert_q_score_to_p_err(*c as f64))
            .collect()
    }
    pub fn convert_p_err_to_q_score(p_err: f64) -> f64 {
        -10.0 * p_err.log10()
    }
    pub fn convert_q_score_to_p_err(q_score: f64) -> f64{
        10.0f64.powf(-0.1 * q_score)
    }
}
