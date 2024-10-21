use crate::dna;
use crate::fx_iterator;
use crate::fx_iterator::FxIterator;
use crate::sequence_read::SequenceRead;
use crate::stats;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};

pub fn fq_check(fq_path: &str, quality_value: u8, ascii_bases: u8) {
    process_fastq_and_check(fq_path, quality_value, ascii_bases)
        .expect("Failed to process fastq file.");
}

fn process_fastq_and_check(fq_path: &str, quality_value: u8, ascii_bases: u8) -> io::Result<()> {
    // let file = File::open(fq_path)?;
    // let reader = BufReader::new(file);

    let mut reads_size = Vec::<usize>::new();
    let mut hashmap_db: HashMap<usize, HashMap<char, usize>> =
        HashMap::from([(0, init_dna_char_hashmap())]); // 0: all

    let mut fq_iter = FxIterator::new(fq_path, "fastq").unwrap();
    while let Some(fx_lines) = fq_iter.next() {
        let mut cur_read = SequenceRead::empty_read();
        cur_read.set_name(fx_lines[0].trim());
        cur_read.set_seq(fx_lines[1].trim());
        cur_read.set_qual(fx_lines[3].trim(), ascii_bases);
        add_sequence_count_into_hashmap(&mut hashmap_db, &cur_read);
        reads_size.push(cur_read.get_seq_len());
    }
    write_fq_stat_to_file("test.txt", &hashmap_db, &reads_size)?;
    Ok(())
}

fn write_fq_stat_to_file(
    path: &str,
    hashmap_db: &HashMap<usize, HashMap<char, usize>>,
    reads_size: &[usize],
) -> io::Result<()> {
    let mut output = File::create(path)?;
    writeln!(output, "{}", get_read_stats_string(reads_size))?;
    writeln!(output, "POS\t#bases\t%A\t%C\t%T\t%G\t%N")?;
    for pos in 0..(hashmap_db.len()) {
        writeln!(output, "{}", get_dna_stats_string(hashmap_db, pos))?;
    }
    Ok(())
}

fn get_dna_stats_string(hashmap_db: &HashMap<usize, HashMap<char, usize>>, pos: usize) -> String {
    let total = get_hasmapdb_value(hashmap_db, pos, 'A')
        + get_hasmapdb_value(hashmap_db, pos, 'T')
        + get_hasmapdb_value(hashmap_db, pos, 'C')
        + get_hasmapdb_value(hashmap_db, pos, 'G')
        + get_hasmapdb_value(hashmap_db, pos, 'N');
    let out = format!(
        "{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
        total,
        100.0 * get_hasmapdb_value(hashmap_db, pos, 'A') as f64 / total as f64,
        100.0 * get_hasmapdb_value(hashmap_db, pos, 'C') as f64 / total as f64,
        100.0 * get_hasmapdb_value(hashmap_db, pos, 'T') as f64 / total as f64,
        100.0 * get_hasmapdb_value(hashmap_db, pos, 'G') as f64 / total as f64,
        100.0 * get_hasmapdb_value(hashmap_db, pos, 'N') as f64 / total as f64
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
    hashmap_db: &HashMap<usize, HashMap<char, usize>>,
    pos: usize,
    key: char,
) -> usize {
    let tmp_hasmap = hashmap_db.get(&pos).expect("should not happened");
    *tmp_hasmap.get(&key).expect("should not happend")
}

fn add_sequence_count_into_hashmap(
    hashmap_db: &mut HashMap<usize, HashMap<char, usize>>,
    read: &SequenceRead,
) {
    for (i, char) in read.get_seq().iter().enumerate() {
        add_dna_char_count_into_pos(0, hashmap_db, char);
        add_dna_char_count_into_pos(i + 1, hashmap_db, char);
    }
}

fn add_dna_char_count_into_pos(
    pos: usize,
    hashmap_db: &mut HashMap<usize, HashMap<char, usize>>,
    char: &dna::DNA,
) {
    match char {
        dna::DNA::A(_) => {
            add_char_count_into_pos(pos, hashmap_db, 'A');
        }
        dna::DNA::T(_) => {
            add_char_count_into_pos(pos, hashmap_db, 'T');
        }
        dna::DNA::C(_) => {
            add_char_count_into_pos(pos, hashmap_db, 'C');
        }
        dna::DNA::G(_) => {
            add_char_count_into_pos(pos, hashmap_db, 'G');
        }
        dna::DNA::N => {
            add_char_count_into_pos(pos, hashmap_db, 'N');
        }
    }
}

fn add_char_count_into_pos(
    pos: usize,
    hashmap_db: &mut HashMap<usize, HashMap<char, usize>>,
    char: char,
) {
    let inner_map = hashmap_db.entry(pos).or_insert_with(init_dna_char_hashmap);
    *inner_map.entry(char).or_insert(0) += 1;
}

fn init_dna_char_hashmap() -> HashMap<char, usize> {
    HashMap::from([('A', 0), ('T', 0), ('C', 0), ('G', 0), ('N', 0)])
}

fn get_read_stats_string(reads_size: &[usize]) -> String {
    let mut sorted_reads_size = reads_size.to_vec();
    sorted_reads_size.sort();
    format!(
        "Length: mean={} ; min={} ; med={} ; max={} ; N50={}",
        stats::average(&sorted_reads_size),
        sorted_reads_size.iter().min().unwrap(),
        stats::median(&sorted_reads_size),
        sorted_reads_size.iter().max().unwrap(),
        stats::n50(&sorted_reads_size)
    )
}

pub fn process_fastq_gz(fq_path: &str, quality_value: u8) -> io::Result<()> {
    let file = File::open(fq_path)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    for line in reader.lines() {
        let line_string = line.expect("not a line");
        println!("print > {}", line_string);
    }
    Ok(())
}
