use crate::io::recordset::{OwnedRecordSet, RecordSetConfig};
use crate::io::fxreader::BatchReader;

use std::thread;

use crossbeam;

pub trait PipelineTask: Send + Sync + Clone + 'static {
    type Partial: Send;
    type Final: Send;
    fn run(&mut self, batch: &mut OwnedRecordSet) -> Result<Self::Partial, std::io::Error>; // need to think about Error type 
    fn merge(&mut self, aggregator: &mut Self::Final, partial: Self::Partial);
}

pub fn run_pipeline<T>(mut reader: BatchReader, config: RecordSetConfig, num_wthread: usize, num_queue: usize, mut task: T, mut final_result: T::Final) -> Option<T::Final>
where T: PipelineTask
{
    let (tx_pool, rx_pool) = crossbeam::channel::bounded::<OwnedRecordSet>(num_queue);
    let (tx_work, rx_work) = crossbeam::channel::bounded::<OwnedRecordSet>(num_queue);
    let (tx_result, rx_result) = crossbeam::channel::unbounded::<T::Partial>();
    
    for _ in 0..num_queue {
        tx_pool.send(OwnedRecordSet::new(config.clone())).unwrap();
    }

    // reader thread
    let reader_handle = thread::spawn(move || {
        while let Ok(mut batch) = rx_pool.recv() {
            match reader.fill_batch(&mut batch) {
                Ok(true) => tx_work.send(batch).unwrap(),
                Ok(false) => break,
                Err(e) => panic!("BatchReader fill_batch() Err: {}", e),
            }
        }
    });
    
    // worker thread
    let mut workers_handle: Vec::<thread::JoinHandle<()>> = Vec::with_capacity(num_wthread);
    for _ in 0..num_wthread {
        let w_rx = rx_work.clone();
        let p_tx = tx_pool.clone();
        let r_tx = tx_result.clone();
        let mut t_task = task.clone();

        workers_handle.push(thread::spawn(move || {
            while let Ok(mut batch) = w_rx.recv() {
                let p_result = t_task.run(&mut batch);
                match p_result {
                    Ok(res) => r_tx.send(res).unwrap(),
                    Err(e) => panic!("PipelineTask run() Error: {}", e),
                }
                let _ = p_tx.send(batch).unwrap();
            }
        }));
    }
    drop(tx_result);

    // merge result from worker threads
    while let Ok(res) = rx_result.recv() {
        task.merge(&mut final_result, res);
    }

    // wait for reader thread finish
    if let Err(e) = reader_handle.join() {
        eprintln!("Reader Thread Error: {:?}", e);
    }

    for wh in workers_handle {
        wh.join().unwrap();
    }

    Some(final_result)
}
