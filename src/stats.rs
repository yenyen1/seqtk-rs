use std::collections::HashMap;

#[derive(Debug)]
pub struct Q2PConverter(HashMap<u8, f64>);
impl Q2PConverter {
    pub fn new(asciibase: u8) -> Self {
        let table = (0..180)
            .map(|i| (i as u8 + asciibase, convert_q_score_to_p_err(i as f64)))
            .collect::<HashMap<_, _>>();
        Q2PConverter(table)
    }
    pub fn get_prob(&self, idx: u8) -> f64 {
        *self.0.get(&idx).expect("Index out of range")
    }
}

pub fn convert_p_err_to_q_score(p_err: f64) -> f64 {
    -10.0 * p_err.log10()
}
pub fn convert_q_score_to_p_err(q_score: f64) -> f64 {
    10.0f64.powf(-q_score / 10.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qp_converter() {
        let table = Q2PConverter::new(0);
        assert_eq!(table.get_prob(0), 1.0);
        assert_eq!(table.get_prob(13), convert_q_score_to_p_err(13.0));
        assert_eq!(table.get_prob(90), convert_q_score_to_p_err(90.0));
    }
}
