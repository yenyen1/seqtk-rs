pub struct SequenceRead<'a> {
    // name: String,
    seq: &'a str,
    qual: Vec<u8>,
}

impl<'a> SequenceRead<'a> {
    /// [Construct SequenceRead]
    pub fn create_read(seq_str: &'a str, qual_str: &str, ascii_base: u8) -> SequenceRead<'a> {
        SequenceRead {
            // name: String::new(),
            seq: seq_str,
            qual: qual_str
                .as_bytes()
                .iter()
                .map(|&byte| byte.wrapping_sub(ascii_base))
                .collect(),
        }
    }
    pub fn get_seq(&self) -> &str {
        self.seq
    }
    // pub fn get_name(&self) -> &str {
    //     self.name.split_whitespace().next().unwrap_or("")
    // }
    pub fn get_q_score_vec(&self) -> &[u8] {
        &self.qual
    }
    pub fn convert_p_err_to_q_score(p_err: f64) -> f64 {
        -10.0 * p_err.log10()
    }
    pub fn convert_q_score_to_p_err(q_score: f64) -> f64 {
        10.0f64.powf(-0.1 * q_score)
    }
}
