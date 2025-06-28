pub fn convert_p_err_to_q_score(p_err: f64) -> f64 {
    -10.0 * p_err.log10()
}
pub fn convert_q_score_to_p_err(q_score: f64) -> f64 {
    10.0f64.powf(-0.1 * q_score)
}
