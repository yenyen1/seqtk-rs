use crate::io_utils::{FaReader, FqReader, FxWriter};
use crate::sub_cli::SampleArgs;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Parses FASTQ/A file and sampling according to the seed and fraction.
/// Outputs the results to [`std::io::stdout()`].
///
/// # Arguments
///
/// Check the arguments by `--help`
///
/// # Errors
///
/// Return an error if the operation cannot be completed.
pub fn subsample_fastx(
    fx_path: &str,
    sparas: &SampleArgs,
    is_fasta: bool,
) -> Result<(), std::io::Error> {
    let rand_seed = sparas.random_seed.unwrap_or(11) as u64;
    let sampling_frac = sparas.sample_fraction.unwrap_or(1.0);
    let mut rng = StdRng::seed_from_u64(rand_seed);

    if is_fasta {
        let fa_iter = FaReader::new(fx_path)?;
        let mut fx_writer = FxWriter::new(is_fasta);
        for record in fa_iter.records() {
            if rng.random::<f64>() <= sampling_frac {
                let read = record.unwrap();
                fx_writer.write(read.id(), read.seq(), None, &[])?;
            }
        }
    } else {
        let fq_iter = FqReader::new(fx_path)?;
        let mut fx_writer = FxWriter::new(is_fasta);
        for record in fq_iter.records() {
            if rng.random::<f64>() <= sampling_frac {
                let read = record.unwrap();
                fx_writer.write(read.id(), read.seq(), read.desc(), read.qual())?;
            }
        }
    }
    Ok(())
}
