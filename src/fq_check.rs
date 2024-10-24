use crate::dna::DNA;
use crate::fx_iterator::FxIterator;
use crate::sequence_read::SequenceRead;
use crate::stats;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, prelude::*};

use strum::IntoEnumIterator;

pub fn fq_check(fq_path: &str, quality_value: u8, ascii_bases: u8) {
    if ascii_bases != 33 && ascii_bases != 64 {
        println!("Warning: Input ascii base {} is not 33 or 64.", ascii_bases);
    }
    process_fastq_and_check(fq_path, quality_value, ascii_bases)
        .expect("Failed to process fastq file.");
}

fn process_fastq_and_check(fq_path: &str, quality_value: u8, ascii_bases: u8) -> io::Result<()> {
    // let mut reads_size = Vec::<usize>::new();
    let mut hashmap_db: HashMap<usize, HashMap<DNA, usize>> =
        HashMap::from([(0, init_dna_char_hashmap())]); // 0: Total
    let mut qual_value_hashset: HashSet<char> = HashSet::new();

    let fq_iter = FxIterator::new(fq_path, "fastq").unwrap();
    // while let Some(fx_lines) = fq_iter.next() {
    for fx_lines in fq_iter {
        let mut cur_read = SequenceRead::empty_read();
        // cur_read.set_name(fx_lines[0].trim());
        cur_read.set_seq(fx_lines[1].trim());
        cur_read.set_qual(fx_lines[3].trim(), ascii_bases);
        add_sequence_count_into_hashmap(&mut hashmap_db, &cur_read);
        qual_value_hashset.extend(cur_read.get_qaul_char());
        // reads_size.push(cur_read.get_seq_len());
        // let cur_qual_u8 = cur_read.get_q_qual_score_u8();
        // println!("qual u8: {:?}", cur_qual_u8);
        // println!("avg qual u8: {}", stats::average(&cur_qual_u8));
        // let cur_p_error_f32 = cur_read.get_p_error_score_f32();
        // let q_err = stats::average(&cur_p_error_f32);
        // println!("p err f32: {:?}", cur_p_error_f32);
        // println!("avg q err: {}", q_err);
        // let p_err_to_q = -10.0 * f32::log10(q_err);
        // println!("p_err_to_q: {}", p_err_to_q);
        // break;
    }
    write_fq_stat_to_file("test.txt", &hashmap_db)?;
    Ok(())
}

fn write_fq_stat_to_file(
    path: &str,
    hashmap_db: &HashMap<usize, HashMap<DNA, usize>>,
) -> io::Result<()> {
    let mut output = File::create(path)?;
    // writeln!(output, "{}", get_read_stats_string(reads_size))?;
    writeln!(output, "POS\t#bases\t%A\t%T\t%C\t%G\t%N")?;
    for pos in 0..(hashmap_db.len()) {
        writeln!(output, "{}", get_dna_stats_string(hashmap_db, pos))?;
    }
    Ok(())
}

fn get_dna_stats_string(hashmap_db: &HashMap<usize, HashMap<DNA, usize>>, pos: usize) -> String {
    let total = get_hasmapdb_value(hashmap_db, pos, DNA::A)
        + get_hasmapdb_value(hashmap_db, pos, DNA::T)
        + get_hasmapdb_value(hashmap_db, pos, DNA::C)
        + get_hasmapdb_value(hashmap_db, pos, DNA::G)
        + get_hasmapdb_value(hashmap_db, pos, DNA::N);
    let out = format!(
        "{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
        total,
        100.0 * get_hasmapdb_value(hashmap_db, pos, DNA::A) as f64 / total as f64,
        100.0 * get_hasmapdb_value(hashmap_db, pos, DNA::C) as f64 / total as f64,
        100.0 * get_hasmapdb_value(hashmap_db, pos, DNA::T) as f64 / total as f64,
        100.0 * get_hasmapdb_value(hashmap_db, pos, DNA::G) as f64 / total as f64,
        100.0 * get_hasmapdb_value(hashmap_db, pos, DNA::N) as f64 / total as f64
    );
    match pos {
        0 => {
            format!("Total\t{}", out)
        }
        _ => {
            format!("{}\t{}", pos, out)
        }
    }
}

fn get_hasmapdb_value(
    hashmap_db: &HashMap<usize, HashMap<DNA, usize>>,
    pos: usize,
    key: DNA,
) -> usize {
    let tmp_hasmap = hashmap_db.get(&pos).expect("should not happened");
    *tmp_hasmap.get(&key).expect("should not happend")
}

fn add_sequence_count_into_hashmap(
    hashmap_db: &mut HashMap<usize, HashMap<DNA, usize>>,
    read: &SequenceRead,
) {
    for (i, char) in read.get_seq().iter().enumerate() {
        add_dna_char_count_into_pos(0, hashmap_db, char);
        add_dna_char_count_into_pos(i + 1, hashmap_db, char);
    }
}

// fn add_dna_char_count_into_pos(
//     pos: usize,
//     hashmap_db: &mut HashMap<usize, HashMap<char, usize>>,
//     char: &dna::DNA,
// ) {
//     match char {
//         dna::DNA::A => {
//             add_char_count_into_pos(pos, hashmap_db, 'A');
//         }
//         dna::DNA::T => {
//             add_char_count_into_pos(pos, hashmap_db, 'T');
//         }
//         dna::DNA::C => {
//             add_char_count_into_pos(pos, hashmap_db, 'C');
//         }
//         dna::DNA::G => {
//             add_char_count_into_pos(pos, hashmap_db, 'G');
//         }
//         dna::DNA::N => {
//             add_char_count_into_pos(pos, hashmap_db, 'N');
//         }
//     }
// }

fn add_dna_char_count_into_pos(
    pos: usize,
    hashmap_db: &mut HashMap<usize, HashMap<DNA, usize>>,
    char: &DNA,
) {
    let inner_map = hashmap_db.entry(pos).or_insert_with(init_dna_char_hashmap);
    *inner_map.entry(char.clone()).or_insert(0) += 1;
}

fn init_dna_char_hashmap() -> HashMap<DNA, usize> {
    let mut dna_hashmap = HashMap::new();
    for c in DNA::iter() {
        dna_hashmap.insert(c, 0);
    }
    dna_hashmap
    // HashMap::from([('A', 0), ('T', 0), ('C', 0), ('G', 0), ('N', 0)])
}

fn get_read_stats_string(reads_size: &[usize]) -> String {
    let mut sorted_reads_size = reads_size.to_vec();
    sorted_reads_size.sort();
    // format!(
    //     "Length: mean={} ; min={} ; med={} ; max={} ; N50={}",
    //     // stats::average(&sorted_reads_size),
    //     sorted_reads_size.iter().min().unwrap(),
    //     stats::median(&sorted_reads_size),
    //     sorted_reads_size.iter().max().unwrap(),
    //     stats::n50(&sorted_reads_size)
    // )
    " ".to_string()
}
