// use rayon::iter::ParallelDrainRange;
// use rayon::prelude::*;
use ndarray::{Array1, Array2};

// use crate::dna::DNA;
use crate::sequence_read::SequenceRead;
// use std::fmt::format;

pub enum QualMap {
    /// [ Qual count map ]
    /// If qual threshold is zero, print out the number of qual count for each position.
    /// If qual threshold is not zero, print out the % of low/high qual for each position.
    DNACount(Array2<u32>), // array2[pos][dna] = count
    QualCount(Array2<u32>),    // array2[pos][qual] = count
    LowQualCount(Array1<u32>), // array1[pos] = low_qual_count

    /// [ Low qual sum ]
    /// Calc avgQ and errQ for each position.
    QualSum(Array2<f64>), // array2[Pos][P/Q] = sum of Qual (P=0 ; Q=1)
}

impl QualMap {
    /// [Initial fqchk map]
    /// Initial DNA count matrix array2[pos][dna] = count
    pub fn init_dna_count_mat(max_read_length: usize) -> QualMap {
        QualMap::DNACount(Array2::<u32>::zeros((max_read_length, 5)))
    }
    /// Initial Qual sum map {ScoreType: {Pos: sum of Qual}}
    pub fn init_qual_sum_mat(max_read_length: usize) -> QualMap {
        QualMap::QualSum(Array2::<f64>::zeros((max_read_length, 2)))
    }
    /// Initial QualMap depend on the qual threshold.
    pub fn init_qual_count_map(
        max_read_length: usize,
        uniq_qual_size: usize,
        qual_threshold: u8,
    ) -> QualMap {
        match qual_threshold {
            0 => {
                let mat = QualMap::init_pos_qual_count(max_read_length, uniq_qual_size);
                QualMap::QualCount(mat)
            }
            _ => {
                let mat = QualMap::init_low_qual_count(max_read_length);
                QualMap::LowQualCount(mat)
            }
        }
    }
    fn init_pos_qual_count(max_read_length: usize, uniq_qual_size: usize) -> Array2<u32> {
        Array2::<u32>::zeros((max_read_length, uniq_qual_size))
    }
    fn init_low_qual_count(max_read_length: usize) -> Array1<u32> {
        Array1::<u32>::zeros(max_read_length)
    }

    /// [ Update data to map ]
    /// Update data with a sequence read
    pub fn update_map_with_read(qual_mat: &mut QualMap, read: &SequenceRead, qual_threshold: u8) {
        match qual_mat {
            QualMap::DNACount(mat) => {
                QualMap::update_dna_count(mat, read);
            }
            QualMap::QualSum(mat) => {
                QualMap::update_qual_sum(mat, read);
            }
            QualMap::QualCount(mat) => {
                QualMap::update_qual_count(mat, read);
            }
            QualMap::LowQualCount(mat) => {
                QualMap::update_low_qual_count(mat, read, qual_threshold);
            }
        }
    }
    fn update_dna_count(dna_count_mat: &mut Array2<u32>, read: &SequenceRead) {
        for (pos, &idx) in read.get_seq_idx().as_slice().iter().enumerate() {
            if let Some(value) = dna_count_mat.get_mut((pos, idx)) {
                *value += 1;
            }
            // dna_count_mat[(pos,idx)] += 1;
        }
    }
    fn update_qual_sum(qual_map: &mut Array2<f64>, read: &SequenceRead) {
        fn add(
            qual_map: &mut Array2<f64>,
            qual_vec: &[f64],
            type_idx: usize, // P=0; Q=1
        ) {
            for (i, &v) in qual_vec.iter().enumerate() {
                if let Some(value) = qual_map.get_mut((type_idx, i)) {
                    *value += v;
                }
            }
        }
        let p_err_vec = read.get_p_err_vec(); // P_error
        add(qual_map, &p_err_vec, 0);
        let q_score_vec = read.get_q_score_vec_f64(); // Q score
        add(qual_map, &q_score_vec, 1);
    }
    fn update_low_qual_count(qual_mat: &mut Array1<u32>, read: &SequenceRead, qual_threshold: u8) {
        for (i, &q) in read.get_q_score_vec().as_slice().iter().enumerate() {
            if q < qual_threshold {
                if let Some(value) = qual_mat.get_mut(i) {
                    *value += 1;
                }
            }
        }
    }
    fn update_qual_count(qual_map: &mut Array2<u32>, read: &SequenceRead) {
        // need to add
    }

    /// [Stats output]
    pub fn get_output(
        dna_count_map: &QualMap,
        qual_sum_map: &QualMap,
        // qual_count_map: &QualMap,
        pos: usize,
    ) -> String {
        let (pos_total, output_dna_stats) =
            QualMap::get_pos_total_and_dna_proportion(dna_count_map, pos);
        let qual_sum_stats = QualMap::get_qual_sum_stats(qual_sum_map, pos, pos_total);
        // let qual_count_stats = QualMap::get_qual_count_stats(qual_count_map,pos, pos_total);
        let mut out = String::new();
        out.push_str(&output_dna_stats);
        out.push_str(&qual_sum_stats);
        // out.push_str(&qual_count_stats);
        out
    }
    pub fn get_qual_count_stats(qual_map: &QualMap, pos: usize, total: f64) -> String {
        match qual_map {
            QualMap::LowQualCount(map) => QualMap::get_low_qual_stats(map, pos, total),
            QualMap::QualCount(_) => String::new(),
            _ => panic!("1 should not happen"),
        }
    }
    /// Get \t{%Low}\t{%High} String from low_qual_count_map
    fn get_low_qual_stats(low_qual_count_mat: &Array1<u32>, pos: usize, total: f64) -> String {
        let low_qual_total = low_qual_count_mat[pos];
        let p_low = 100.0 * low_qual_total as f64 / total;
        format!("\t{:.1}\t{:.1}", p_low, 100.0 - p_low)
    }
    /// Get Total\t{#Base}\t{%A}\t{%C}\t{%G}\t{%T}\t{%N} String from dna_count_mat
    pub fn get_total_and_dna_proportion(dna_count_mat: &QualMap) -> (f64, String) {
        match dna_count_mat {
            QualMap::DNACount(mat) => {
                let total_dna_counts = mat.sum_axis(ndarray::Axis(0)); // Sum of count for A/C/G/T/N (across rows/pos)
                let total = total_dna_counts.sum() as f64;
                let prop_str: String = total_dna_counts
                    .iter()
                    .map(|&c| format!("\t{:.1}", 100.0 * c as f64 / total))
                    .collect();
                (total, format!("Total\t{:.0}{}", total, prop_str))
            }
            _ => panic!("1 should not happen"),
        }
    }
    /// Get {Pos}\t{#Base}\t{%A}\t{%C}\t{%G}\t{%T}\t{%N} String from dna_count_mat
    pub fn get_pos_total_and_dna_proportion(dna_count_mat: &QualMap, pos: usize) -> (f64, String) {
        match dna_count_mat {
            QualMap::DNACount(mat) => {
                let pos_total: u32 = mat.slice(ndarray::s![pos, ..]).iter().sum(); // Sum of counts at position
                let total_f64: f64 = pos_total as f64;
                let prop_str: String = mat
                    .slice(ndarray::s![pos, ..])
                    .iter()
                    .map(|&c| format!("\t{:.1}", 100.0 * c as f64 / total_f64))
                    .collect();
                (total_f64, format!("{}\t{}{}", pos + 1, pos_total, prop_str))
            }
            _ => panic!("2 should not happen"),
        }
    }
    /// Get {avgQ}\t{Qerr} String from qual_sum_map
    pub fn get_qual_sum_stats(qual_sum_map: &QualMap, pos: usize, total: f64) -> String {
        match qual_sum_map {
            QualMap::QualSum(mat) => {
                // P=0; Q=1
                let q_avg = mat.get((pos, 1)).unwrap_or(&0.0) / total;
                let p_avg = mat.get((pos, 0)).unwrap_or(&0.0) / total;
                let q_err = SequenceRead::convert_p_err_to_q_score(p_avg);
                format!("\t{:.1}\t{:.1}", q_avg, q_err)
            }
            _ => {
                panic!("3 should not happen")
            }
        }
    }
}
