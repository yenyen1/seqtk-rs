use rayon::slice::ParallelSliceMut;

use crate::io::fxreader::{BatchReader, FxReader};
use crate::io::parallel::{run_pipeline, PipelineTask};
use crate::io::recordset::{OwnedRecordSet, RecordSetConfig};

use std::path::Path;
use std::process::ExitCode;

#[derive(Debug, Clone)]
pub struct BatchSize {
    read_count: usize,
    base_count: u64,
    min_length: u32,
    max_length: u32,
    median_length: f64,
    average: f64,
    n50: u32,
    lengths: Vec<u32>,
}
impl BatchSize {
    pub fn new() -> Self {
        // the same as record_size_limit of RecordSetConfig
        Self::with_capacity(16 * 1024)
    }
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            read_count: 0,
            base_count: 0,
            min_length: 0,
            max_length: 0,
            median_length: 0.0,
            average: 0.0,
            n50: 0,
            lengths: Vec::with_capacity(capacity),
        }
    }
    pub fn add(&mut self, length: u32) {
        // self.read_count += 1;
        self.base_count += length as u64;
        self.lengths.push(length);
    }
    pub fn merge(&mut self, size: BatchSize) {
        // self.read_count += size.read_count;
        self.base_count += size.base_count;
        self.lengths.extend_from_slice(&size.lengths);
    }
    pub fn calculate_statistics(&mut self) {
        self.lengths.par_sort_unstable();
        self.read_count = self.lengths.len();

        if self.read_count > 0 {
            self.min_length = self.lengths[0];
            self.max_length = self.lengths[self.read_count - 1];
            self.average = self.base_count as f64 / self.read_count as f64;

            let mid = self.read_count / 2;
            self.median_length = match self.read_count % 2 {
                1 => self.lengths[mid] as f64,
                _ => (self.lengths[mid - 1] + self.lengths[mid]) as f64 / 2.0,
            };

            let half: u64 = self.base_count / 2;
            let mut acc: u64 = 0;
            for &cur in self.lengths.iter().rev() {
                acc += cur as u64;
                if acc >= half {
                    self.n50 = cur;
                    break;
                }
            }
        };
    }
}
impl std::fmt::Display for BatchSize {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // #seq, #bases, avg_size, min_size, med_size, max_size, N50
        writeln!(
            f,
            "{:?}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.read_count,
            self.base_count,
            self.average,
            self.min_length,
            self.median_length,
            self.max_length,
            self.n50
        )
    }
}
impl Default for BatchSize {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone)]
pub struct SizeTask;

impl PipelineTask for SizeTask {
    type Partial = BatchSize;
    type Final = BatchSize;

    fn run(&mut self, batch: &mut OwnedRecordSet) -> BatchSize {
        let mut batch_size = BatchSize::new();
        for record in batch.iter_seq() {
            batch_size.add(record.seq.len() as u32);
        }
        batch_size
    }
    fn merge(&mut self, aggregator: &mut BatchSize, partial: BatchSize) {
        aggregator.merge(partial);
    }
}

pub fn run<P: AsRef<Path>>(path: P) -> ExitCode {
    let batch_reader = match FxReader::new(path) {
        Ok(r) => BatchReader::new(r),
        Err(e) => {
            eprintln!("[ERROR] {}", e);
            return ExitCode::from(2);
        }
    };
    let recordset_config = RecordSetConfig::fasta_default();

    let result = run_pipeline(
        batch_reader,
        recordset_config,
        2,
        4,
        SizeTask,
        BatchSize::with_capacity(10 * 1024 * 1024),
    );
    match result {
        Ok(mut r) => {
            r.calculate_statistics();
            println!("{}", r);
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("{}", e);
            ExitCode::from(2)
        }
    }
}
