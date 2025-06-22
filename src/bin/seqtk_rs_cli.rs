use clap::Parser;
use seqtk_rs::{fq_check, seq, sub_cli, subsample};
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
            let sample_paras = subsample::SampleParas::new(
                sample.random_seed.unwrap_or(11),
                sample.sample_fraction,
            );
            if let Some(in_fq) = &sample.in_fq {
                subsample::subsample_fastx(in_fq, &sample_paras, false)?;
            }
            if let Some(in_fa) = &sample.in_fa {
                subsample::subsample_fastx(in_fa, &sample_paras, true)?;
            }
        }

        sub_cli::Commands::Seq(seq) => {
            sub_cli::valiation_seq_args(seq)?;
            let filter_rule = seq::FilterParas::from(seq);
            let mask_paras = seq::MaskParas::from(seq);
            let out_paras = seq::OutArgs::from(seq);
            if let Some(in_fq) = &seq.in_fq {
                seq::parse_fastx(in_fq, &filter_rule, &mask_paras, &out_paras, false)?;
            }
            if let Some(in_fa) = &seq.in_fa {
                seq::parse_fastx(in_fa, &filter_rule, &mask_paras, &out_paras, true)?;
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
