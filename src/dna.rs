// ASCII to ACGTN
// A|a: 0 ; C|c:1 ; G|g:2 ; T|t:3 ; other: 4
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
pub fn get_dna_idx_from_u8(s: u8) -> usize {
    ASCII_TO_ACGTN_IDX[s as usize]
}

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
    /// ( a: 7; c: 8; g: 9; t:10 ) ; ( ryswkm:11 ) ; ( bdhv:12 ) ; ( n:13 ) ; others: 14
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
    fn get_iupac_index_from_u8(s: u8) -> usize {
        Self::ASCII_TO_IUPAC_NC_IDX[s as usize]
    }
    /// Input IUPAC index and output 1: CG, 2: CG(reverse complemet), 3: others
    fn get_all_cp_index(cur_b: usize, next_b: usize) -> usize {
        if (cur_b == 1 || cur_b == 8) && (next_b == 2 || next_b == 9) {
            1
        } else if (cur_b == 2 || cur_b == 9) && (next_b == 1 || next_b == 8) {
            2
        } else {
            3
        }
    }
    /// Input IUPAC index and output 1: CG, 2: CG(reverse complemet), 3: others
    fn get_unmasked_cp_index(cur_b: usize, next_b: usize) -> usize {
        if cur_b == 1 && next_b == 2 {
            1
        } else if cur_b == 2 && next_b == 1 {
            2
        } else {
            3
        }
    }
    pub fn get_all_result(v: &[usize; 23]) -> [usize; 9] {
        // #A, #C, #G, #T, #2, #3, #4, #CpG, #GC
        let mut out: [usize; 9] = [0; 9];
        out[0] = v[0] + v[7];
        out[1] = v[1] + v[8];
        out[2] = v[2] + v[9];
        out[3] = v[3] + v[10];
        out[4] = v[4] + v[11];
        out[5] = v[5] + v[12];
        out[6] = v[6] + v[13];
        out[7] = v[15];
        out[8] = v[16];
        out
    }
    pub fn get_unmasked_result(v: &[usize; 23]) -> [usize; 9] {
        // #A, #C, #G, #T, #2, #3, #4, #CG, #GC
        let mut out: [usize; 9] = [0; 9];
        out[0] = v[0];
        out[1] = v[1];
        out[2] = v[2];
        out[3] = v[3];
        out[4] = v[4];
        out[5] = v[5];
        out[6] = v[6];
        out[7] = v[15];
        out[8] = v[16];
        out
    }
    pub fn count_all_nc(count: &mut [usize; 23], seq: &[u8], start: usize, end: usize) {
        // count: 0-14 'A'-'X'; 15-17 CpG, CpG-rc, others
        // let mut count: [usize; 23] = [0; 23];
        let (mut cur_b, mut next_b): (usize, usize) = (0, 0);
        seq[start..end].windows(2).for_each(|b| {
            cur_b = Self::get_iupac_index_from_u8(b[0]);
            next_b = Self::get_iupac_index_from_u8(b[1]);
            count[cur_b] += 1;
            count[Self::get_all_cp_index(cur_b, next_b) + 14] += 1;
        });
        cur_b = Self::get_iupac_index_from_u8(seq[end - 1]);
        count[cur_b] += 1;
    }
    pub fn count_unmasked_nc(count: &mut [usize; 23], seq: &[u8], start: usize, end: usize) {
        // let mut count: [usize; 23] = [0; 23];
        let mut cur_b: usize = 0;
        seq[start..end].windows(2).for_each(|b| {
            cur_b = Self::get_iupac_index_from_u8(b[0]);
            if cur_b < 7 {
                count[cur_b] += 1;
                let next_b = Self::get_iupac_index_from_u8(b[1]);
                if next_b < 7 {
                    count[Self::get_unmasked_cp_index(cur_b, next_b) + 14] += 1;
                }
            }
        });
        let cur_b = Self::get_iupac_index_from_u8(seq[end - 1]);
        if cur_b < 7 {
            count[cur_b] += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;

    #[test]
    fn test_get_cp_index() {
        fn get_cp(cur_b: u8, next_b: u8, excluded_masked: bool) -> usize {
            if excluded_masked {
                SeqComp::get_unmasked_cp_index(
                    SeqComp::get_iupac_index_from_u8(cur_b),
                    SeqComp::get_iupac_index_from_u8(next_b),
                )
            } else {
                SeqComp::get_all_cp_index(
                    SeqComp::get_iupac_index_from_u8(cur_b),
                    SeqComp::get_iupac_index_from_u8(next_b),
                )
            }
        }
        assert_eq!(get_cp(b'C', b'G', false), 1, "err1");
        assert_eq!(get_cp(b'c', b'G', false), 1, "err2");
        assert_eq!(get_cp(b'G', b'C', false), 2, "err3");
        assert_eq!(get_cp(b'g', b'c', false), 2, "err4");
        assert_eq!(get_cp(b'g', b'g', false), 3, "err5");
        assert_eq!(get_cp(b'n', b'g', false), 3, "err6");
        assert_eq!(get_cp(b'a', b'T', false), 3, "err7");

        assert_eq!(get_cp(b'C', b'G', true), 1, "err8");
        assert_eq!(get_cp(b'c', b'g', true), 3, "err9");
        assert_eq!(get_cp(b'G', b'C', true), 2, "err10");
        assert_eq!(get_cp(b'g', b'c', true), 3, "err11");
        assert_eq!(get_cp(b'g', b'G', true), 3, "err12");
        assert_eq!(get_cp(b'C', b'C', true), 3, "err13");
        assert_eq!(get_cp(b'A', b'T', true), 3, "err14");
        assert_eq!(get_cp(b'N', b'T', true), 3, "err15");
    }

    #[test]
    fn test_count_nucleotides() {
        fn get_count(seq: &[u8], start: usize, end: usize, all: bool) -> [usize; 9] {
            let mut count: [usize; 23] = [0; 23];
            if all {
                SeqComp::count_all_nc(&mut count, seq, start, end);
                SeqComp::get_all_result(&count)
            } else {
                SeqComp::count_unmasked_nc(&mut count, seq, start, end);
                SeqComp::get_unmasked_result(&count)
            }
        }
        let record = Record::with_attrs(
            "ID1",
            None,
            b"CGTACCGCGACGATCGATcgACTTGCGTNNnNAGTCRYSWKMCCCATCgCTTryswkmBDhV",
            b"qewr4gjiroeggfryremb[trdgqewr4gjiroeggfryremb[trdgryremb[trdge",
        );
        let count = get_count(record.seq(), 0, 25, true);
        assert_eq!(count, [5, 8, 7, 5, 0, 0, 0, 6, 1], "err1");
        let count = get_count(record.seq(), 25, 58, true);
        assert_eq!(count, [2, 7, 3, 5, 12, 0, 4, 2, 1], "err2");
        let count = get_count(record.seq(), 58, 62, true);
        assert_eq!(count, [0, 0, 0, 0, 0, 4, 0, 0, 0], "err3");
        let count = get_count(record.seq(), 0, record.seq().len(), true);
        assert_eq!(count, [7, 15, 10, 10, 12, 4, 4, 8, 3], "err4");

        let count = get_count(record.seq(), 0, 25, false);
        assert_eq!(count, [5, 7, 6, 5, 0, 0, 0, 5, 1], "err5");
        let count = get_count(record.seq(), 25, 58, false);
        assert_eq!(count, [2, 7, 2, 5, 6, 0, 3, 1, 0], "err6");
        let count = get_count(record.seq(), 58, 62, false);
        assert_eq!(count, [0, 0, 0, 0, 0, 3, 0, 0, 0], "err7");
        let count = get_count(record.seq(), 0, record.seq().len(), false);
        assert_eq!(count, [7, 14, 8, 10, 6, 3, 3, 6, 2], "err8");
    }
}

// / IUPAC nucleotide code: "ACMGRSVTWYHKDBN"
// / https://www.bioinformatics.org/sms/iupac.html
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
