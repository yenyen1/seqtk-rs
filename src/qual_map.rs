use crate::sequence_read::SequenceRead;

use ndarray::{Array1, Array2, Axis};
use rayon::prelude::*;
use std::collections::HashMap;

pub enum QualMap {
    /// [ Qual count map ]
    /// If qual threshold is zero, print out the number of qual count for each position.
    QualCount(Array2<f64>), // array2[pos][qual] = count
    /// If qual threshold is not zero, print out the % of low/high qual for each position.
    LowQualCount(Array1<f64>), // array1[pos] = low_qual_count
}

impl QualMap {
    /// Initial QualMap depend on the qual threshold.
    pub fn init_qual_count_arr(
        max_read_len: usize,
        qual_size: usize,
        qual_threshold: u8,
    ) -> QualMap {
        match qual_threshold {
            0 => {
                let mat = Array2::<f64>::zeros((max_read_len, qual_size));
                QualMap::QualCount(mat)
            }
            _ => {
                let mat = Array1::<f64>::zeros(max_read_len);
                QualMap::LowQualCount(mat)
            }
        }
    }

    /// [ Update data to map ]
    /// Update data with a sequence read
    pub fn update_qual_count_mat(
        &mut self,
        read: &SequenceRead,
        qual_idx_map: &HashMap<u8, usize>,
        qual_threshold: u8,
    ) {
        match self {
            QualMap::QualCount(mat) => {
                Self::update_qual_count(mat, read, qual_idx_map);
            }
            QualMap::LowQualCount(mat) => {
                Self::update_low_qual_count(mat, read, qual_threshold);
            }
        }
    }
    fn update_low_qual_count(inner_mat: &mut Array1<f64>, read: &SequenceRead, qual_threshold: u8) {
        for (i, &q) in read.get_q_score_vec().iter().enumerate() {
            if q < qual_threshold {
                if let Some(value) = inner_mat.get_mut(i) {
                    *value += 1.0;
                }
            }
        }
    }
    fn update_qual_count(
        inner_mat: &mut Array2<f64>,
        read: &SequenceRead,
        qual_idx_map: &HashMap<u8, usize>,
    ) {
        for (i, q) in read.get_q_score_vec().iter().enumerate() {
            let qual_idx = qual_idx_map.get(q).expect("qual_map-066 should not happen");
            if let Some(value) = inner_mat.get_mut((i, *qual_idx)) {
                *value += 1.0;
            }
        }
    }

    /// [Stats output]
    pub fn get_qual_count_stats(&self, pos: usize, total: f64) -> String {
        match self {
            QualMap::LowQualCount(mat) => Self::get_low_qual_stats(mat, pos, total),
            QualMap::QualCount(mat) => Self::get_qual_stats(mat, pos, total),
        }
    }
    pub fn get_total_qual_count_stats(&self, total: f64) -> String {
        match self {
            QualMap::QualCount(mat) => Self::get_total_qual_stats(mat, total),
            QualMap::LowQualCount(mat) => Self::get_total_low_qual_stats(mat, total),
        }
    }
    /// Get \t{%Q1}\t{%Q2}...
    fn get_qual_stats(qual_mat: &Array2<f64>, pos: usize, total: f64) -> String {
        let mut out = String::new();
        for &q in qual_mat.slice(ndarray::s![pos, ..]) {
            out.push_str(&format!("\t{:.1}", 100.0 * q / total));
        }
        out
    }
    /// Get \t{%Low}\t{%High} String from low_qual_count_mat
    fn get_low_qual_stats(low_qual_mat: &Array1<f64>, pos: usize, total: f64) -> String {
        let p_low = 100.0 * low_qual_mat.get(pos).unwrap_or(&0.0) / total;
        format!("\t{:.1}\t{:.1}", p_low, 100.0 - p_low)
    }
    /// Get \t{%Q1}\t{%Q2}...
    fn get_total_qual_stats(qual_mat: &Array2<f64>, total: f64) -> String {
        let total_qual: Vec<f64> = qual_mat
            .axis_iter(Axis(1)) // Sum across rows
            .into_par_iter()
            .map(|col| col.sum())
            .collect();
        let mut out = String::new();
        for &q in total_qual.iter() {
            out.push_str(&format!("\t{:.1}", 100.0 * q / total));
        }
        out
    }
    /// Get \t{%Low}\t{%High} String from low_qual_count_mat
    fn get_total_low_qual_stats(low_qual_mat: &Array1<f64>, total: f64) -> String {
        let p_low = 100.0 * low_qual_mat.sum() / total;
        format!("\t{:.1}\t{:.1}", p_low, 100.0 - p_low)
    }
}
