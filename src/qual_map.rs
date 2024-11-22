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
    QualCount(Array2<u32>), // array2[pos][qual] = count
    LowQualCount(Array1<u32>), // array1[pos] = low_qual_count
}

impl QualMap {
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
            QualMap::QualCount(mat) => {
                QualMap::update_qual_count(mat, read);
            }
            QualMap::LowQualCount(mat) => {
                QualMap::update_low_qual_count(mat, read, qual_threshold);
            }
        }
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
}
