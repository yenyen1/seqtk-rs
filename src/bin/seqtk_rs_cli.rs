use clap::Parser;
use seqtk_rs::{fqchk, nc_comp, seq, size, sub_cli, subsample, trim};

/// Here 1
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = sub_cli::Cli::parse();

    match &cli.command {
        sub_cli::Commands::Fqchk(fqchk) => {
            let qthreshold = fqchk.quality_value.unwrap_or(0);
            let ascii = fqchk.ascii_base.unwrap_or(33);
            if qthreshold == 0 {
                fqchk::get_result_wo_qthreshold(&fqchk.in_fq, ascii as usize)?;
            } else {
                fqchk::get_result_with_qthreshold(
                    &fqchk.in_fq,
                    qthreshold + ascii,
                    ascii as usize,
                )?;
            }
        }

        sub_cli::Commands::Sample(sample) => {
            if let Some(fq) = &sample.in_fq {
                subsample::subsample_fastx(fq, sample, false)?;
            }
            if let Some(fa) = &sample.in_fa {
                subsample::subsample_fastx(fa, sample, true)?;
            }
        }

        sub_cli::Commands::Size(size) => {
            if let Some(fq) = &size.in_fq {
                size::calc_fq_size(fq)?;
            }
            if let Some(fa) = &size.in_fa {
                size::calc_fa_size(fa)?;
            }
        }

        sub_cli::Commands::Qctrim(trim) => {
            let q_threshold = trim.q_thershold.unwrap_or(13);
            let asciibase = trim.ascii_base.unwrap_or(33);
            let minlen = trim.min_length.unwrap_or(30);
            trim::trimfq(&trim.in_fq, q_threshold + asciibase, minlen)?;
        }
        sub_cli::Commands::Comp(comp) => {
            if let Some(fq) = &comp.in_fq {
                match &comp.in_bed {
                    Some(bed) => nc_comp::calc_fq_comp_with_bed(fq, bed, comp.exclude_masked)?,
                    None => nc_comp::calc_fq_comp_wo_bed(fq, comp.exclude_masked)?,
                }
            }
            if let Some(fa) = &comp.in_fa {
                match &comp.in_bed {
                    Some(bed) => nc_comp::calc_fa_comp_with_bed(fa, bed, comp.exclude_masked)?,
                    None => nc_comp::calc_fa_comp_wo_bed(fa, comp.exclude_masked)?,
                }
            }
        }

        sub_cli::Commands::Seq(seq) => {
            sub_cli::valiation_seq_args(seq)?;
            if let Some(path) = &seq.in_fq {
                seq::parse_fastq(path, seq)?;
            }
            if let Some(path) = &seq.in_fa {
                seq::parse_fasta(path, seq)?;
            }
        }
    }
    Ok(())
}
