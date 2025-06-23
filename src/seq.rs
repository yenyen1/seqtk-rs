use crate::bed::{self, BedMap, BedPos};
use crate::dna;
use crate::io::{FaReader, FqReader, FxWriter};
use crate::record::RecordType;
use crate::sub_cli::SeqArgs;

struct FilterParas {
    mini_seq_length: usize,
    drop_ambigous_seq: bool,
    output_odd_reads: bool,
    output_even_reads: bool,
}
impl FilterParas {
    fn from(seq: &SeqArgs) -> Self {
        FilterParas {
            mini_seq_length: seq.mini_seq_length.unwrap_or(0),
            drop_ambigous_seq: seq.drop_ambigous_seq,
            output_odd_reads: seq.output_even,
            output_even_reads: seq.output_odd,
        }
    }
}
struct MaskParas {
    mask_char: Option<char>,
    uppercases: bool,
    lowercases_to_char: bool,
    q_low: u8,
    q_high: u8,
    mask_regions: Option<String>,
    mask_complement_region: bool,
}
impl MaskParas {
    fn from(seq: &SeqArgs) -> Self {
        let ascii_bases = seq.ascii_bases.unwrap_or(33);
        MaskParas {
            mask_char: seq.mask_char,
            uppercases: seq.uppercases,
            lowercases_to_char: seq.lowercases_to_char,
            q_low: (seq.q_low.unwrap_or(0) + ascii_bases),
            q_high: (seq.q_high.unwrap_or(255 - ascii_bases) + ascii_bases),
            mask_regions: seq.mask_regions.clone(),
            mask_complement_region: seq.mask_complement_region,
        }
    }
}
pub struct OutArgs {
    output_qual_shift: u8,
    fake_fastq_quality: Option<char>,
    output_fasta: bool,
    reverse_complement: bool,
    both_complement: bool,
    trim_header: bool,
    line_len: Option<usize>,
}
impl OutArgs {
    fn from(seq: &SeqArgs) -> Self {
        let ascii_bases = seq.ascii_bases.unwrap_or(33);
        let out_qual_shift = if seq.output_qual_33 {
            ascii_bases - 33
        } else {
            0
        };
        let output_fasta =
            if !seq.output_fasta && seq.in_fa.is_some() && seq.fake_fastq_quality.is_none() {
                true
            } else {
                seq.output_fasta
            };
        OutArgs {
            output_qual_shift: out_qual_shift,
            fake_fastq_quality: seq.fake_fastq_quality,
            output_fasta,
            reverse_complement: seq.reverse_complement,
            both_complement: seq.both_complement,
            trim_header: seq.trim_header,
            line_len: seq.line_len,
        }
    }
}

pub fn parse_fasta(path: &str, seq: &SeqArgs) -> Result<(), std::io::Error> {
    let fparas = FilterParas::from(seq);
    let mparas = MaskParas::from(seq);
    let oparas = OutArgs::from(seq);
    let bed_map = match &mparas.mask_regions {
        Some(bed_path) => bed::get_bed_map(bed_path)?,
        None => BedMap::new(),
    };
    let fa_iter = FaReader::new(path)?;
    let mut fx_writer = FxWriter::new(oparas.output_fasta);
    for (i, record) in fa_iter.records().enumerate() {
        let read = record.unwrap();
        let is_pass = is_pass(i + 1, &read, &fparas);
        if is_pass {
            modify_and_print_read(&mut fx_writer, &read, &mparas, &oparas, &bed_map, true)?;
        }
    }
    Ok(())
}
pub fn parse_fastq(path: &str, seq: &SeqArgs) -> Result<(), std::io::Error> {
    let fparas = FilterParas::from(seq);
    let mparas = MaskParas::from(seq);
    let oparas = OutArgs::from(seq);
    let bed_map = match &mparas.mask_regions {
        Some(bed_path) => bed::get_bed_map(bed_path)?,
        None => BedMap::new(),
    };

    let fq_iter = FqReader::new(path)?;
    let mut fx_writer = FxWriter::new(oparas.output_fasta);
    for (i, record) in fq_iter.records().enumerate() {
        let read = record.unwrap();
        let is_pass = is_pass(i + 1, &read, &fparas);
        if is_pass {
            modify_and_print_read(&mut fx_writer, &read, &mparas, &oparas, &bed_map, false)?;
        }
    }

    Ok(())
}
fn modify_and_print_read(
    fx_writer: &mut FxWriter,
    read: &dyn RecordType,
    mask_paras: &MaskParas,
    out_paras: &OutArgs,
    bed_map: &BedMap,
    is_fasta: bool,
) -> Result<(), std::io::Error> {
    let mut seq = modify_seq(read, mask_paras, bed_map, is_fasta);
    let desc = if out_paras.trim_header {
        None
    } else {
        read.desc()
    };
    let mut qual = modify_qual(read, out_paras, is_fasta);

    if let Some(line_len) = out_paras.line_len {
        add_newlines(&mut seq, line_len);
        add_newlines(&mut qual, line_len);
    }

    if out_paras.reverse_complement {
        revcomp(&mut seq, &mut qual);
        fx_writer.write(read.id(), &seq, desc, &qual)?;
    } else if out_paras.both_complement {
        fx_writer.write(read.id(), &seq, desc, &qual)?;
        revcomp(&mut seq, &mut qual);
        fx_writer.write(read.id(), &seq, desc, &qual)?;
    } else {
        fx_writer.write(read.id(), &seq, desc, &qual)?;
    }
    Ok(())
}
fn add_newlines(data: &mut Vec<u8>, line_len: usize) {
    let mut i = line_len;
    while i < data.len() {
        data.insert(i, b'\n'); // Insert '\n' after i
        i += line_len + 1;
    }
}
fn revcomp(seq: &mut [u8], qual: &mut [u8]) {
    dna::revcomp(seq);
    qual.reverse();
}
fn modify_qual(read: &dyn RecordType, oparas: &OutArgs, is_fasta: bool) -> Vec<u8> {
    match oparas.fake_fastq_quality {
        Some(fake_qual) => {
            vec![fake_qual as u8; read.seq().len()]
        }
        None => {
            if is_fasta {
                Vec::new()
            } else if oparas.output_qual_shift == 0 {
                read.qual().to_vec()
            } else {
                read.qual()
                    .iter()
                    .map(|q| q - oparas.output_qual_shift)
                    .collect()
            }
        }
    }
}
fn modify_seq(
    read: &dyn RecordType,
    mask_paras: &MaskParas,
    bed_map: &BedMap,
    is_fasta: bool,
) -> Vec<u8> {
    let default_bed_pos = vec![BedPos(usize::MAX, 0)];
    let bed_pos = bed_map.get(read.id()).unwrap_or(&default_bed_pos);
    let mut seq = if mask_paras.uppercases {
        read.seq().to_ascii_uppercase() // convert all to uppercases
    } else {
        read.seq().to_vec()
    };

    match mask_paras.mask_char {
        Some(c) => {
            let c_u8 = c as u8;
            let mut cur_idx: usize = 0;
            seq.iter_mut().enumerate().for_each(|(i, ch)| {
                if mask_paras.lowercases_to_char && ch.is_ascii_lowercase() {
                    *ch = c_u8;
                } else {
                    let (is_overlap, c_idx) = bed::is_overlapping(i, &bed_pos[cur_idx..]);
                    cur_idx += c_idx;
                    if (is_overlap && !mask_paras.mask_complement_region)
                        || (!is_overlap && mask_paras.mask_complement_region)
                    {
                        *ch = c_u8;
                    }
                }
            })
        }
        None => {
            let mut cur_idx: usize = 0;
            seq.iter_mut().enumerate().for_each(|(i, ch)| {
                let (is_overlap, c_idx) = bed::is_overlapping(i, &bed_pos[cur_idx..]);
                cur_idx += c_idx;
                if (is_overlap && !mask_paras.mask_complement_region)
                    || (!is_overlap && mask_paras.mask_complement_region)
                {
                    *ch = ch.to_ascii_lowercase();
                }
            })
        }
    }
    if !is_fasta {
        match mask_paras.mask_char {
            Some(c) => {
                let c_u8 = c as u8;
                read.qual()
                    .iter()
                    .zip(seq.iter_mut())
                    .for_each(|(&qual, ch)| {
                        if qual < mask_paras.q_low || mask_paras.q_high < qual {
                            *ch = c_u8; // mask bases by q_low and q_high
                        }
                    })
            }
            None => {
                read.qual()
                    .iter()
                    .zip(seq.iter_mut())
                    .for_each(|(&qual, ch)| {
                        if qual < mask_paras.q_low || mask_paras.q_high < qual {
                            *ch = ch.to_ascii_lowercase(); // mask bases by q_low and q_high
                        }
                    })
            }
        }
    }
    seq
}
fn is_pass(i: usize, read: &dyn RecordType, fparas: &FilterParas) -> bool {
    if fparas.output_even_reads && i % 2 == 0 {
        return false;
    }
    if fparas.output_odd_reads && i % 2 == 1 {
        return false;
    }
    if fparas.drop_ambigous_seq && read.seq().iter().any(|&c| dna::get_dna_idx_from_u8(c) > 3) {
        return false;
    }
    if read.seq().len().le(&fparas.mini_seq_length) {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;

    #[test]
    fn test_modify_qual() {
        fn init_oparas() -> OutArgs {
            OutArgs {
                output_qual_shift: 0,
                fake_fastq_quality: None,
                output_fasta: false,
                reverse_complement: false,
                both_complement: false,
                trim_header: false,
                line_len: None,
            }
        }
        let record = Record::with_attrs("SEQ_ID_1", None, b"ATCGATcgACTTG", b"gfryremb[trdg");
        let mut oparas = init_oparas();

        // [01] without modify
        let out_qual = modify_qual(&record, &oparas, false);
        assert_eq!(&out_qual, b"gfryremb[trdg");

        // [02] check output_qual_shift
        oparas.output_qual_shift = 10;
        let out_qual = modify_qual(&record, &oparas, false);
        assert_eq!(&out_qual, b"]\\hoh[cXQjhZ]");

        // [02] check output_qual_shift
        oparas = init_oparas();
        oparas.fake_fastq_quality = Some('T');
        let out_qual = modify_qual(&record, &oparas, false);
        assert_eq!(&out_qual, b"TTTTTTTTTTTTT");
    }

    #[test]
    fn test_modify_seq() {
        fn init_mparas() -> MaskParas {
            MaskParas {
                mask_char: None,
                uppercases: false,
                lowercases_to_char: false,
                q_low: 0,
                q_high: 255,
                mask_regions: None,
                mask_complement_region: false,
            }
        }
        let record = Record::with_attrs("SEQ_ID_1", None, b"ATCGATcgACTTG", b"!(*AAAABbbaaz");
        let mut bed_map: BedMap = BedMap::new();
        let mut mparas = init_mparas();

        // [01] without modify
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"ATCGATcgACTTG", "[err01]");

        // [02] check mask_char + lowrcases_to_char
        mparas.mask_char = Some('N');
        mparas.lowercases_to_char = true;
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"ATCGATNNACTTG", "[err02]");

        // [03] check uppercases
        mparas = init_mparas();
        mparas.uppercases = true;
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"ATCGATCGACTTG", "[err03]");

        // [04] check q_low and q_high
        mparas = init_mparas();
        mparas.q_low = 20 + 33;
        mparas.q_high = 255;
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"atcGATcgACTTG", "[err04]");
        mparas.q_low = 0;
        mparas.q_high = 85 + 33;
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"ATCGATcgACTTg", "[err05-1]");

        // [05] check q_low and q_high + mask_char and lowercases_to_char
        mparas = init_mparas();
        mparas.mask_char = Some('N');
        mparas.q_low = 20 + 33;
        mparas.q_high = 85 + 33;
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"NNNGATcgACTTN", "[err05-2]");
        mparas.lowercases_to_char = true;
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"NNNGATNNACTTN", "[err05-3]");

        // [06] check q_low and q_high + uppercases
        mparas = init_mparas();
        mparas.uppercases = true;
        mparas.q_low = 20 + 33;
        mparas.q_high = 85 + 33;
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"atcGATCGACTTg", "[err06]");

        // [07] check mask by bed
        let bed_path = Some("test.bed".to_string());
        bed_map.add("SEQ_ID_1".to_string(), BedPos(5, 10));
        mparas = init_mparas();
        mparas.mask_regions = bed_path.clone();
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"ATCGAtcgacTTG", "[err07]");

        // [08] checl mask by complement region
        mparas.mask_complement_region = true;
        let out_seq = modify_seq(&record, &mparas, &bed_map, false);
        assert_eq!(&out_seq, b"atcgaTcgACttg", "[err08]");
    }

    #[test]
    fn test_is_filtered() {
        fn init_fparas() -> FilterParas {
            FilterParas {
                mini_seq_length: 0,
                drop_ambigous_seq: false,
                output_odd_reads: false,
                output_even_reads: false,
            }
        }
        let mut fparas = init_fparas();
        let record = Record::with_attrs("@SEQ_ID_1", None, b"ATCGATCGACTTG", b"!!<AAAABbbaab");

        // [01] check mini_seq_length
        fparas.mini_seq_length = 10;
        assert_eq!(is_pass(0, &record, &fparas), true);
        fparas.mini_seq_length = 50;
        assert_eq!(is_pass(0, &record, &fparas), false);

        // [02] check drop ambiguous seq
        fparas = init_fparas();
        fparas.drop_ambigous_seq = true;
        assert_eq!(is_pass(0, &record, &fparas), true);
        let record2 = Record::with_attrs("@SEQ_ID_2", None, b"ANCGATCGACTTG", b"!!<AAAABbbaab");
        assert_eq!(is_pass(0, &record2, &fparas), false);

        // [03] output odd read
        fparas.output_odd_reads = true;
        assert_eq!(is_pass(0, &record, &fparas), true);
        assert_eq!(is_pass(1, &record, &fparas), false);

        // [04] output even read
        fparas = init_fparas();
        fparas.output_even_reads = true;
        assert_eq!(is_pass(3, &record, &fparas), true);
        assert_eq!(is_pass(4, &record, &fparas), false);
    }
    #[test]
    fn test_add_newlines() {
        let mut seq = b"aaaaabbbbbcccccdddddeeeeefffff".to_vec();
        add_newlines(&mut seq, 5);
        assert_eq!(seq, b"aaaaa\nbbbbb\nccccc\nddddd\neeeee\nfffff");
    }
    #[test]
    fn test_revcomp() {
        let mut seq = b"AATTCCGG".to_vec();
        let mut qual = b"<<((vv++".to_vec();

        revcomp(&mut seq, &mut qual);

        assert_eq!(seq, b"CCGGAATT");
        assert_eq!(qual, b"++vv((<<");
    }
}
