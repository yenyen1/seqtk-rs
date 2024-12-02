use clap::{Args, Parser, Subcommand};
use seqtk_rs::fq_check;
use std::time::Instant;

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// fastq QC summary
    Fqchk(FqchkArgs),
    /// Seq
    Seq(SeqArgs),
}

#[derive(Args)]
struct FqchkArgs {
    /// input fastq path
    in_fq: String,
    /// output tsv path
    #[arg(short, long)]
    out: Option<String>,

    #[arg(short, long)]
    /// quality value [default: 0]
    quality_value: Option<u8>,

    #[arg(short, long)]
    /// ascii value [default: 33]
    ascii_base: Option<u8>,
}

#[derive(Args)]
struct SeqArgs {
    #[arg(long)]
    /// input fastq path or fasta path
    in_fq: Option<String>,
    /// out path
    #[arg(long)]
    out: Option<String>,

    #[arg(long)]
    /// mask bases that quality lower than q_low [default: 0]
    q_low: Option<u8>,
    /// mask bases that quality higher than q_high [default: 255]
    #[arg(long)]
    q_high: Option<u8>,
    #[arg(long)]
    /// mask bases converted to CHAR [default: convert to lowercase]
    mask_char: Option<char>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let cli = Cli::parse();

    match &cli.command {
        Commands::Fqchk(fqchk) => {
            fq_check::fq_check(
                &fqchk.in_fq,
                &fqchk.out.clone().unwrap_or("test".to_string()),
                fqchk.quality_value.unwrap_or(0),
                fqchk.ascii_base.unwrap_or(33),
            )?;
        }
        Commands::Seq(seq) => {}
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
