use crate::io::{FaReader, FqReader, FxWriter};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

pub struct SampleParas {
    random_seed: u64,
    sample_fraction: Option<f64>,
}
impl SampleParas {
    pub fn new(random_seed: u64, sample_fraction: Option<f64>) -> Self {
        SampleParas {
            random_seed,
            sample_fraction,
        }
    }
}

pub fn subsample_fastx(
    fx_path: &str,
    sample_paras: &SampleParas,
    is_fasta: bool,
) -> Result<(), std::io::Error> {
    let sampling_frac = sample_paras.sample_fraction.unwrap_or(1.0);
    let mut rng = StdRng::seed_from_u64(sample_paras.random_seed);

    if is_fasta {
        let fa_iter = FaReader::new(fx_path)?;
        let mut fx_writer = FxWriter::new(is_fasta);
        for record in fa_iter.records() {
            if rng.random::<f64>() > sampling_frac {
                continue;
            }
            let read = record.unwrap();
            fx_writer.write(read.id(), read.seq(), None, &[])?;
        }
    } else {
        let fq_iter = FqReader::new(fx_path)?;
        let mut fx_writer = FxWriter::new(is_fasta);
        for record in fq_iter.records() {
            if rng.random::<f64>() > sampling_frac {
                continue;
            }
            let read = record.unwrap();
            fx_writer.write(read.id(), read.seq(), read.desc(), read.qual())?;
        }
    }
    Ok(())
}
