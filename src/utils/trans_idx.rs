const fn build_dna_table() -> [u8; 256] {
    let mut table = [4u8; 256];

    let mut i = 0;
    while i < 256 {
        let b = i as u8;
        table[i] = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        };
        i += 1;
    }
    table
}

// get TABLE in compile time
pub const ASCII_TO_ACGTN_IDX: [u8; 256] = build_dna_table();
