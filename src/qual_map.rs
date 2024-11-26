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
                for (i, q) in read.get_q_score_vec().iter().enumerate() {
                    let qual_idx = qual_idx_map.get(q).expect("qual_map-066 should not happen");
                    mat[(i, *qual_idx)] += 1.0;
                }
            }
            QualMap::LowQualCount(mat) => {
                for (i, &q) in read.get_q_score_vec().iter().enumerate() {
                    if q < qual_threshold {
                        mat[i] += 1.0;
                    }
                }
            }
        }
    }
    /// [Stats output]
    pub fn get_qual_count_stats(&self, pos: usize, total: f64) -> String {
        match self {
            QualMap::LowQualCount(mat) => {
                // Get \t{%Low}\t{%High} String from low_qual_count_mat
                let p_low = 100.0 * mat[pos] / total;
                format!("\t{:.1}\t{:.1}", p_low, 100.0 - p_low)
            },
            QualMap::QualCount(mat) => {
                // Get \t{%Q1}\t{%Q2}...
                let mut out = String::new();
                for &q in mat.slice(ndarray::s![pos, ..]) {
                    out.push_str(&format!("\t{:.1}", 100.0 * q / total));
                }
                out
            },
        }
    }
    pub fn get_total_qual_count_stats(&self, total: f64) -> String {
        match self {
            QualMap::QualCount(mat) => {
                // Get \t{%Q1}\t{%Q2}...
                let total_qual: Vec<f64> = mat
                    .axis_iter(Axis(1)) // Sum across rows
                    .into_par_iter()
                    .map(|col| col.sum())
                    .collect();
                let mut out = String::new();
                for &q in total_qual.iter() {
                    out.push_str(&format!("\t{:.1}", 100.0 * q / total));
                }
                out
            },
            QualMap::LowQualCount(mat) => {
                // Get \t{%Low}\t{%High} String from low_qual_count_mat
                let p_low = 100.0 * mat.sum() / total;
                format!("\t{:.1}\t{:.1}", p_low, 100.0 - p_low)
            },
        }
    }
}
