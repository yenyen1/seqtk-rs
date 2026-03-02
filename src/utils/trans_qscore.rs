
pub struct QScoreTable {
    p_err_table: [f32; 256],
}
impl QScoreTable {
    pub fn new(ascii_base: u8) -> Self {
        let mut p_err_table = [1.0f32; 256];
        for val in ascii_base..=255u8 {
            let q_score = val.saturating_sub(ascii_base) as f32;
            p_err_table[val as usize] = 10f32.powf(-q_score/10.0);
        }
        Self {
            p_err_table,
        }
    }
    #[inline(always)]
    pub fn get_p_err(&self, ascii_char: u8) -> f32{
        self.p_err_table[ascii_char as usize]
    }
}