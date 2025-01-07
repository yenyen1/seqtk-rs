use crate::dna;
use crate::utils::{self, FxWriter, RecordType};
use bio::io::fastq::Record;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::collections::HashMap;

pub struct FilterParas {
    // fastq
    mini_seq_length: usize,       // fq(V), fa(need to add)
    drop_ambigous_seq: bool,      // fq(V), fa(need to add)
    output_odd_reads: bool,       // fq, fa
    output_even_reads: bool,      // fq, fa
    random_seed: u64,             // fq
    sample_fraction: Option<f64>, // fq
}
impl FilterParas {
    pub fn new(
        mini_seq_length: usize,
        drop_ambigous_seq: bool,
        output_odd_reads: bool,
        output_even_reads: bool,
        random_seed: u64,
        sample_fraction: Option<f64>,
    ) -> Self {
        assert!(
            !(output_even_reads && output_odd_reads),
            "Should not use --output-even-reads and --output-odd-reads options together."
        );
        FilterParas {
            mini_seq_length,
            drop_ambigous_seq,
            output_odd_reads,
            output_even_reads,
            random_seed,
            sample_fraction,
        }
    }
}
pub struct MaskParas<'a> {
    mask_char: Option<char>,
    uppercases: bool,
    lowercases_to_char: bool,
    q_low: u8,                        // fq
    q_high: u8,                       // fq
    mask_regions: &'a Option<String>, // fa
    mask_complement_region: bool,     // fa
}
impl<'a> MaskParas<'a> {
    pub fn new(
        mask_char: Option<char>,
        uppercases: bool,
        lowercases_to_char: bool,
        q_low: u8,
        q_high: u8,
        mask_regions: &'a Option<String>,
        mask_complement_region: bool,
    ) -> Self {
        assert!(
            !(mask_complement_region && mask_regions.is_none()),
            "--mask-complment-region should effective with --mask-regions."
        );
        assert!(
            !(lowercases_to_char && mask_char.is_none()),
            "--lowercases-to-char should effective with --mask-char."
        );
        MaskParas {
            mask_char,
            uppercases,
            lowercases_to_char,
            q_low,
            q_high,
            mask_regions,
            mask_complement_region,
        }
    }
}
pub struct OutArgs {
    output_qual_shift: u8, // fq Q
    fake_fastq_quality: Option<char>,
    output_fasta: bool,       //
    reverse_complement: bool, // fa | fq
    both_complement: bool,    // fa | fq
    trim_header: bool,
    line_len: Option<usize>,
}
impl OutArgs {
    pub fn new(
        output_qual_shift: u8, // fq Q
        fake_fastq_quality: Option<char>,
        output_fasta: bool,
        reverse_complement: bool,
        both_complement: bool,
        trim_header: bool,
        line_len: Option<usize>,
    ) -> Self {
        assert!(
            !(output_fasta && (output_qual_shift != 0 || fake_fastq_quality.is_some())),
            "Should not use --output-fasta with --output-qual-shift or --fake-fastq-quality."
        );
        assert!(
            !(output_qual_shift != 0 && fake_fastq_quality.is_some()),
            "Should not use --output-qual-shift and --fake-fastq-quality together."
        );
        assert!(
            !(reverse_complement && both_complement),
            "Should not use --reverse-complement and --both-complement together."
        );
        OutArgs {
            output_qual_shift,
            fake_fastq_quality,
            output_fasta,
            reverse_complement,
            both_complement,
            trim_header,
            line_len,
        }
    }
}

pub fn parse_fastx(
    fx_path: &str,
    filter_rule: &FilterParas,
    mask_paras: &MaskParas,
    out_paras: &OutArgs,
    is_fasta: bool,
) -> Result<(), std::io::Error> {
    let bed_map = match mask_paras.mask_regions {
        Some(bed_path) => utils::get_bed_map(bed_path)?,
        None => HashMap::new(),
    };
    let mut fx_writer = FxWriter::new(out_paras.output_fasta);
    match filter_rule.sample_fraction {
        Some(frac) => {
            let mut rng = StdRng::seed_from_u64(filter_rule.random_seed);
            if is_fasta {
                let mut fa_iter = utils::new_fa_iterator(fx_path)?;
            } else {
                let mut fq_iter = utils::new_fq_iterator(fx_path)?;
                for (i, record) in fq_iter.records().enumerate() {
                    let read = record.unwrap();
                    let filter = filter_read(i + 1, &read, filter_rule);
                    if filter && rng.gen::<f64>() <= frac {
                        modify_and_print_read(
                            &mut fx_writer,
                            &read,
                            mask_paras,
                            out_paras,
                            &bed_map,
                        )?;
                    }
                }
            }
        }
        None => {
            let fx_iter = utils::new_fq_iterator(fx_path)?;
            if is_fasta {
                let mut fa_iter = utils::new_fa_iterator(fx_path)?;
            } else {
                for (i, record) in fx_iter.records().enumerate() {
                    let read = record.unwrap();
                    let filtered_read = filter_read(i + 1, &read, filter_rule);
                    if filtered_read {
                        modify_and_print_read(
                            &mut fx_writer,
                            &read,
                            mask_paras,
                            out_paras,
                            &bed_map,
                        )?;
                    }
                }
            }
        }
    }
    Ok(())
}
fn modify_and_print_read(
    fx_writer: &mut FxWriter,
    read: &dyn RecordType,
    mask_paras: &MaskParas,
    out_paras: &OutArgs,
    bed_map: &HashMap<String, Vec<[usize; 2]>>,
) -> Result<(), std::io::Error> {
    let mut seq = modify_seq(read, mask_paras, bed_map);
    let desc = if out_paras.trim_header {
        None
    } else {
        read.desc()
    };
    let mut qual = modify_qual(read, out_paras);

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
fn modify_qual(read: &dyn RecordType, out_paras: &OutArgs) -> Vec<u8> {
    match out_paras.fake_fastq_quality {
        Some(fake_qual) => {
            vec![fake_qual as u8; read.seq().len()]
        }
        None => {
            if out_paras.output_qual_shift == 0 {
                read.qual().to_vec()
            } else {
                read.qual()
                    .iter()
                    .map(|q| q - out_paras.output_qual_shift)
                    .collect()
            }
        }
    }
}
fn modify_seq(
    read: &dyn RecordType,
    mask_paras: &MaskParas,
    bed_map: &HashMap<String, Vec<[usize; 2]>>,
) -> Vec<u8> {
    fn is_overlapping(pos: usize, bed_pos: &Vec<[usize; 2]>) -> bool {
        bed_pos
            .iter()
            .any(|&range| range[0] <= pos && pos < range[1])
    }
    let default_bed_pos = vec![[usize::MAX, 0]];
    let bed_pos = bed_map.get(read.id()).unwrap_or(&default_bed_pos);
    let mut seq = if mask_paras.uppercases {
        read.seq().to_ascii_uppercase() // convert all to uppercases
    } else {
        read.seq().to_vec()
    };
    match mask_paras.mask_char {
        Some(c) => {
            let c_u8 = c as u8;
            if mask_paras.lowercases_to_char {
                for (i, (&qual, ch)) in read.qual().iter().zip(seq.iter_mut()).enumerate() {
                    if qual < mask_paras.q_low || mask_paras.q_high < qual {
                        *ch = c_u8; // mask bases by q_low and q_high
                    } else if ch.is_ascii_lowercase() {
                        *ch = c_u8; // convert lowercase to char
                    } else {
                        let is_overlap = is_overlapping(i, bed_pos);
                        if (mask_paras.mask_complement_region && !is_overlap)
                            || (!mask_paras.mask_complement_region && is_overlap)
                        {
                            *ch = c_u8; // mask bases by bed
                        }
                    }
                }
            } else {
                for (i, (&qual, ch)) in read.qual().iter().zip(seq.iter_mut()).enumerate() {
                    if qual < mask_paras.q_low || mask_paras.q_high < qual {
                        *ch = c_u8; // mask bases by q_low and q_high
                    } else {
                        let is_overlap = is_overlapping(i, bed_pos);
                        if (mask_paras.mask_complement_region && !is_overlap)
                            || (!mask_paras.mask_complement_region && is_overlap)
                        {
                            *ch = c_u8; // mask bases by bed
                        }
                    }
                }
            }
        }
        None => {
            for (i, (&qual, ch)) in read.qual().iter().zip(seq.iter_mut()).enumerate() {
                if qual < mask_paras.q_low || mask_paras.q_high < qual {
                    *ch = ch.to_ascii_lowercase(); // mask bases by q_low and q_high
                } else {
                    let is_overlap = is_overlapping(i, bed_pos);
                    if (mask_paras.mask_complement_region && !is_overlap)
                        || (!mask_paras.mask_complement_region && is_overlap)
                    {
                        *ch = ch.to_ascii_lowercase(); // mask bases by bed
                    }
                }
            }
        }
    }
    seq
}
fn filter_read(i: usize, read: &dyn RecordType, filter_rule: &FilterParas) -> bool {
    let even_or_odd = if !(filter_rule.output_even_reads && filter_rule.output_odd_reads) {
        if filter_rule.output_even_reads {
            i % 2 == 0
        } else if filter_rule.output_odd_reads {
            i % 2 == 1
        } else {
            true
        }
    } else {
        false
    };
    let seq_wo_n: bool = if filter_rule.drop_ambigous_seq {
        !read.seq().iter().any(|&c| dna::get_dna_idx_from_u8(c) > 3) // not containing char other than ATCGatcg
    } else {
        true
    };
    let write_read = read.seq().len().gt(&filter_rule.mini_seq_length) && // seq len large than min_len 
        seq_wo_n && // seq not containing N
        even_or_odd;
    write_read
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;

    #[test]
    fn test_modify_qual() {
        let record = Record::with_attrs("SEQ_ID_1", None, b"ATCGATcgACTTG", b"gfryremb[trdg");
        // [01] no modify
        let oparas = OutArgs::new(0, None, false, false, false, false, None);
        let out_qual = modify_qual(&record, &oparas);
        assert_eq!(&out_qual, b"gfryremb[trdg");
        // [02] check output_qual_shift
        let oparas = OutArgs::new(10, None, false, false, false, false, None);
        let out_qual = modify_qual(&record, &oparas);
        assert_eq!(&out_qual, b"]\\hoh[cXQjhZ]");
        // [02] check output_qual_shift
        let oparas = OutArgs::new(0, Some('T'), false, false, false, false, None);
        let out_qual = modify_qual(&record, &oparas);
        assert_eq!(&out_qual, b"TTTTTTTTTTTTT");
    }
    #[test]
    fn test_modify_seq() {
        let record = Record::with_attrs("SEQ_ID_1", None, b"ATCGATcgACTTG", b"!(*AAAABbbaaz");
        let mut bed_map: HashMap<String, Vec<[usize; 2]>> = HashMap::new();

        // [01] no modify
        let mparas = MaskParas::new(None, false, false, 0, 255, &None, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"ATCGATcgACTTG", "[err01]");
        // [02] check mask_char + lowrcases_to_char
        let mparas = MaskParas::new(Some('N'), false, true, 0, 255, &None, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"ATCGATNNACTTG", "[err02]");
        // [03] check uppercases
        let mparas = MaskParas::new(None, true, false, 0, 255, &None, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"ATCGATCGACTTG", "[err03]");
        // [04] check q_low and q_high
        let mparas = MaskParas::new(None, false, false, 20 + 33, 255, &None, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"atcGATcgACTTG", "[err04]");
        let mparas = MaskParas::new(None, false, false, 0, 85 + 33, &None, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"ATCGATcgACTTg", "[err05-1]");
        // [05] check q_low and q_high + mask_char and lowercases_to_char
        let mparas = MaskParas::new(Some('N'), false, false, 20 + 33, 85 + 33, &None, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"NNNGATcgACTTN", "[err05-2]");
        let mparas = MaskParas::new(Some('N'), false, true, 20 + 33, 85 + 33, &None, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"NNNGATNNACTTN", "[err05-3]");
        // [06] check q_low and q_high + uppercases
        let mparas = MaskParas::new(None, true, false, 20 + 33, 85 + 33, &None, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"atcGATCGACTTg", "[err06]");
        // [07] check mask by bed
        let bed_path = Some("test.bed".to_string());
        bed_map.insert("SEQ_ID_1".to_string(), vec![[5, 10]]);

        let mparas = MaskParas::new(None, false, false, 0, 255, &bed_path, false);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"ATCGAtcgacTTG", "[err07]");
        // [08] checl mask by complement region
        let mparas = MaskParas::new(None, false, false, 0, 255, &bed_path, true);
        let out_seq = modify_seq(&record, &mparas, &bed_map);
        assert_eq!(&out_seq, b"atcgaTcgACttg", "[err08]");
    }
    #[test]
    fn test_filter_read() {
        let record = Record::with_attrs("@SEQ_ID_1", None, b"ATCGATCGACTTG", b"!!<AAAABbbaab");

        // [01] check mini_seq_length
        let fparas = FilterParas::new(10, false, false, false, 0, None);
        assert_eq!(filter_read(0, &record, &fparas), true);
        let fparas = FilterParas::new(50, false, false, false, 0, None);
        assert_eq!(filter_read(0, &record, &fparas), false);

        // [02] check drop anbugous seq
        let fparas = FilterParas::new(0, true, false, false, 0, None);
        assert_eq!(filter_read(0, &record, &fparas), true);
        let record2 = Record::with_attrs("@SEQ_ID_2", None, b"ANCGATCGACTTG", b"!!<AAAABbbaab");
        let fparas = FilterParas::new(0, true, false, false, 0, None);
        assert_eq!(filter_read(0, &record2, &fparas), false);

        // [03] output odd read
        let fparas = FilterParas::new(0, true, true, false, 0, None);
        assert_eq!(filter_read(1, &record, &fparas), true);
        let fparas = FilterParas::new(0, true, true, false, 0, None);
        assert_eq!(filter_read(0, &record, &fparas), false);

        // [04] output even read
        let fparas = FilterParas::new(0, true, false, true, 0, None);
        assert_eq!(filter_read(2, &record, &fparas), true);
        let fparas = FilterParas::new(0, true, false, true, 0, None);
        assert_eq!(filter_read(3, &record, &fparas), false);
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
