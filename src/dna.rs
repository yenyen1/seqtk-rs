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
    for c in seq.iter_mut() {
        *c = complement(c);
    }
    seq.reverse();
}
