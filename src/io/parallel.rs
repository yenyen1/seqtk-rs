use crate::io::fxerr::FxError;
use crate::io::fxreader::BatchReader;
use crate::io::recordset::{OwnedRecordSet, RecordSetConfig};

use std::thread;

use crossbeam;

/// Parallel Task for worker threads
pub trait PipelineTask: Send + Sync + Clone + 'static {
    type Partial: Send;
    type Final: Send;
    /// The worker function invoked in parallel to process individual RecordSet batches.
    /// This function performs the core computation on a data subset.
    /// It is expected to be infallible (no errors should occur).
    fn run(&mut self, batch: &mut OwnedRecordSet) -> Self::Partial;

    /// Merges a partial result into the global aggregator to produce the final result.
    /// This is invoked as soon as an individual parallel task completes,
    /// without waiting for other pending tasks to finish.
    fn merge(&mut self, aggregator: &mut Self::Final, partial: Self::Partial);
}

pub fn run_pipeline<T>(
    mut reader: BatchReader,
    config: RecordSetConfig,
    num_wthread: usize,
    num_queue: usize,
    mut task: T,
    mut final_result: T::Final,
) -> Result<T::Final, FxError>
// need to check is there any Error return?
where
    T: PipelineTask,
{
    let (tx_pool, rx_pool) = crossbeam::channel::bounded::<OwnedRecordSet>(num_queue);
    let (tx_work, rx_work) = crossbeam::channel::bounded::<OwnedRecordSet>(num_queue);
    let (tx_result, rx_result) = crossbeam::channel::unbounded::<T::Partial>();

    for _ in 0..num_queue {
        tx_pool.send(OwnedRecordSet::new(config.clone())).unwrap();
    }

    // reader thread
    let reader_handle = thread::spawn(move || -> Result<(), FxError> {
        while let Ok(mut batch) = rx_pool.recv() {
            match reader.fill_batch(&mut batch) {
                Ok(true) => tx_work.send(batch).unwrap(),
                Ok(false) => break,
                // Error return from `seq-io` crate (BufferLimit or FormatError when parsing records)
                Err(e) => {
                    return Err(e);
                }
            }
        }
        Ok(())
    });

    // worker thread
    let mut workers_handle: Vec<thread::JoinHandle<()>> = Vec::with_capacity(num_wthread);
    for _ in 0..num_wthread {
        let w_rx = rx_work.clone();
        let p_tx = tx_pool.clone();
        let r_tx = tx_result.clone();
        let mut t_task = task.clone();

        workers_handle.push(thread::spawn(move || {
            while let Ok(mut batch) = w_rx.recv() {
                let p_result = t_task.run(&mut batch);
                r_tx.send(p_result).unwrap(); // need to recheck Error scenario
                let _ = p_tx.send(batch); // need to recheck Error
            }
        }));
    }
    drop(tx_result);

    // merge result from worker threads
    while let Ok(res) = rx_result.recv() {
        task.merge(&mut final_result, res);
    }

    // wait for reader thread finish
    match reader_handle.join().unwrap() {
        Ok(_) => eprintln!("[INFO] Finished reading File."),
        Err(e) => eprintln!(
            "[ERROR] BatchReader fill_batch() terminated unexpectedly. {}",
            e
        ),
    }

    for wh in workers_handle {
        wh.join().unwrap();
    }

    Ok(final_result)
}
