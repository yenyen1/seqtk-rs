use clap::Parser;
use seqtk_rs::{fq_check, nc_comp, seq, size, sub_cli, subsample};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let cli = sub_cli::Cli::parse();
    let out = "test".to_string();

    match &cli.command {
        sub_cli::Commands::Fqchk(fqchk) => {
            fq_check::fq_check(
                &fqchk.in_fq,
                &fqchk.out.clone().unwrap_or(out),
                fqchk.quality_value.unwrap_or(0),
                fqchk.ascii_base.unwrap_or(33),
            )?;
        }

        sub_cli::Commands::Sample(sample) => {
            if let Some(in_fq) = &sample.in_fq {
                subsample::subsample_fastx(in_fq, sample, false)?;
            }
            if let Some(in_fa) = &sample.in_fa {
                subsample::subsample_fastx(in_fa, sample, true)?;
            }
        }

        sub_cli::Commands::Size(size) => {
            if let Some(in_fq) = &size.in_fq {
                size::calc_fq_size(in_fq)?;
            }
            if let Some(in_fa) = &size.in_fa {
                size::calc_fa_size(in_fa)?;
            }
        }

        sub_cli::Commands::Comp(size) => {
            if let Some(in_fq) = &size.in_fq {
                nc_comp::calc_fq_comp(in_fq)?;
            }
            if let Some(in_fa) = &size.in_fa {
                nc_comp::calc_fa_comp(in_fa)?;
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
    eprintln!(
        "Process time: {:?}",
        format_duration(start.elapsed().as_secs())
    );
    Ok(())
}

fn format_duration(duration: u64) -> String {
    let seconds = duration;
    let minutes = seconds / 60;
    let hours = minutes / 60;
    let days = hours / 24;
    let formatted_time = format!(
        "{}-{:02}:{:02}:{:02}",
        days,
        hours % 24,
        minutes % 60,
        seconds % 60
    );
    formatted_time
}
