use clap::Parser;
use seqtk_rs::{fq_check, seq, sub_cli};
use std::time::Instant;

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    pub command: sub_cli::Commands,
}

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
        sub_cli::Commands::Seq(seq) => {
            let filter_rule = seq::FilterParas::new(
                seq.mini_seq_length.unwrap_or(0),
                seq.drop_ambigous_seq,
                seq.output_even_reads,
                seq.output_odd_reads,
                seq.random_seed.unwrap_or(11),
                seq.sample_fraction,
            );
            let ascii_bases = seq.quality_shift.unwrap_or(33);
            let mask_paras = seq::MaskParas::new(
                seq.q_low.unwrap_or(0) + ascii_bases,
                seq.q_high.unwrap_or(222) + ascii_bases,
                seq.mask_char,
                &seq.mask_regions,
                seq.mask_complement_region,
                seq.uppercases,
                seq.lowercases_to_char,
            );
            seq::parse_seq(
                &seq.in_fx,
                &seq.out.clone().unwrap_or(out),
                &filter_rule,
                &mask_paras,
            )?;
        }
    }
    println!(
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
