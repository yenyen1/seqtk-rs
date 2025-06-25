/// ASCII to ACGTN
/// A|a: 0 ; C|c:1 ; G|g:2 ; T|t:3 ; other: 4
const ASCII_TO_ACGTN_IDX: [usize; 256] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];
/// complemet for IUPAC nucleotide code
/// https://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
const COMP_TRANS: [u8; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, b'T', b'V', b'G', b'H', 69, 70,
    b'C', b'D', 73, 74, b'M', 76, b'K', b'N', 79, 80, 81, b'Y', b'S', b'A', b'A', b'B', b'W', 88,
    b'R', 90, 91, 92, 93, 94, 95, 96, b't', b'v', b'G', b'h', 101, 102, b'c', b'd', 105, 106, b'm',
    108, b'k', b'n', 111, 112, 113, b'y', b's', b'a', b'a', b'b', b'w', 120, b'r', 122, 123, 124,
    125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
    163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181,
    182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200,
    201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
    220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238,
    239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
];
pub const ACGTN_TAB: &str = "\tA\tC\tG\tT\tN";
/// IUPAC nucleotide code: "ACMGRSVTWYHKDBN"
/// https://www.bioinformatics.org/sms/iupac.html
// const IUPAC: [u8; 256] = [
//      0,  1,  2,  3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
// 	16, 17, 18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
// 	32, 33, 34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
// 	48, 49, 50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
//     64, b'A', b'B', b'C', b'D', 69, 70, b'G', b'H', 73, 74, b'K', 76, b'M', b'N', 79, // 64 - 79
//     80, 81, b'R', b'S', b'T', b'U', b'V', b'W', 88, b'Y', 90, 91, 92, 93, 94, 95, // 80 -
//     96, b'a', b'b', b'c', b'd', 101, 102, b'g', b'h', 105, 106, b'k', 108, b'm', b'n', 111,
//     112, 113, b'r', b's', b't', b'u', b'v', b'w', 120, b'y', 122, 123, 124, 125, 126, 127,
//     128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
//     144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
//     160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
//     176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
//     192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
//     208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
//     224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
//     240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
// ];
pub fn get_dna_idx_from_u8(s: u8) -> usize {
    ASCII_TO_ACGTN_IDX[s as usize]
}

pub fn complement(s: &u8) -> u8 {
    COMP_TRANS[*s as usize]
}

pub fn revcomp(seq: &mut [u8]) {
    seq.iter_mut().for_each(|c| *c = complement(c));
    seq.reverse();
}

pub struct SeqComp;
impl SeqComp {
    /// ASCII to IUPAC nc comsidering masked bases
    /// ( A: 0; C: 1; G: 2; T: 3 ) ; ( RYSWKM: 4 ) ; ( BDHV: 5 ) ; ( N: 6 )
    /// ( a: 7; c: 8; g: 9; t:10 ) ; ( ryswkm:11 ) ; ( bdhv:12 ) ; ( n:13 )
    /// others: 14
    const ASCII_TO_IUPAC_NC_IDX: [usize; 256] = [
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, //
        0, 5, 1, 5, 14, 14, 2, 5, 14, 14, 4, 14, 4, 6, 14, 14, 14, 4, 4, 3, 14, 5, 4, 14, 4, 14,
        14, 14, 14, 14, 14, 14, // A-
        7, 12, 8, 12, 14, 14, 9, 12, 14, 14, 11, 14, 11, 13, 14, 14, 14, 11, 11, 10, 14, 12, 11,
        14, 11, 14, 14, 14, 14, 14, 14, 14, // a-
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
    ];
    /// Purine (A|G) = 0; Pyrimidine (C|T) = 1
    const IUPAC_TO_PURINE_PYRIMIDINE: [usize; 15] = [
        0, 1, 0, 1, 2, 2, 2, //
        0, 1, 0, 1, 2, 2, 2, 2,
    ];
    fn get_iupac_index_from_u8(s: u8) -> usize {
        Self::ASCII_TO_IUPAC_NC_IDX[s as usize]
    }
    /// Input IUPAC index and output 1: transversion, 2: transition , 3: others
    fn get_tx_index(cur_b: usize, next_b: usize) -> usize {
        let cur_b = Self::IUPAC_TO_PURINE_PYRIMIDINE[cur_b];
        let next_b = Self::IUPAC_TO_PURINE_PYRIMIDINE[next_b];
        if cur_b == 2 || next_b == 2 {
            3
        } else if cur_b == next_b {
            2
        } else {
            1
        }
    }
    /// Input IUPAC index and output 1: CG, 2: CG(reverse complemet), 3: CG(masked), 4: CG(masked rc), 5: others
    fn get_cp_index(cur_b: usize, next_b: usize) -> usize {
        if cur_b == 1 && next_b == 2 {
            1
        } else if cur_b == 2 && next_b == 1 {
            2
        } else if cur_b == (1 | 8) && next_b == (2 | 9) {
            3
        } else if cur_b == (2 | 9) && next_b == (1 | 8) {
            4
        } else {
            5
        }
    }
    fn get_comp_output(v: &[usize; 23], include_masked: bool) -> [usize; 11] {
        // #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts
        let mut out: [usize; 11] = [0; 11];
        if include_masked {
            out[0] = v[0] + v[7];
            out[1] = v[1] + v[8];
            out[2] = v[2] + v[9];
            out[3] = v[3] + v[10];
            out[4] = v[4] + v[11];
            out[5] = v[5] + v[12];
            out[6] = v[6] + v[13];
            out[7] = v[18] + v[19] + v[20] + v[21];
            out[10] = v[18] + v[20];
        } else {
            out[0] = v[0];
            out[1] = v[1];
            out[2] = v[2];
            out[3] = v[3];
            out[4] = v[4];
            out[5] = v[5];
            out[6] = v[6];
            out[7] = v[18] + v[19];
            out[10] = v[18];
        }
        out[8] = v[15];
        out[9] = v[16];
        out
    }
    pub fn count_all_nucleotides(seq: &[u8], start: usize, end: usize) -> [usize; 11] {
        // count: 0-14 'A'-'X'; 15-17 tv, ts, other; 18-22 CpG, CpG-rc, masked-CpG, masked-CpG-rc, others
        let mut count: [usize; 23] = [0; 23];
        seq[start..end].windows(2).for_each(|b| {
            let cur_b = Self::get_iupac_index_from_u8(b[0]);
            let next_b = Self::get_iupac_index_from_u8(b[1]);
            count[cur_b] += 1;
            count[Self::get_tx_index(cur_b, next_b) + 14] += 1;
            count[Self::get_cp_index(cur_b, next_b) + 17] += 1;
        });
        let cur_b = Self::get_iupac_index_from_u8(seq[end - 1]);
        count[cur_b] += 1;
        Self::get_comp_output(&count, true)
    }
    pub fn count_unmasked_nucleotides(seq: &[u8], start: usize, end: usize) -> [usize; 11] {
        let mut count: [usize; 23] = [0; 23];
        seq[start..end].windows(2).for_each(|b| {
            let cur_b = Self::get_iupac_index_from_u8(b[0]);
            if cur_b < 7 {
                count[cur_b] += 1;
                let next_b = Self::get_iupac_index_from_u8(b[1]);
                if next_b < 7 {
                    count[Self::get_tx_index(cur_b, next_b) + 14] += 1;
                    count[Self::get_cp_index(cur_b, next_b) + 17] += 1;
                }
            }
        });
        let cur_b = Self::get_iupac_index_from_u8(seq[end - 1]);
        if cur_b < 7 {
            count[cur_b] += 1;
        }
        Self::get_comp_output(&count, false)
    }
}
