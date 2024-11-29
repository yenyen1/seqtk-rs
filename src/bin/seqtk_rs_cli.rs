use clap::{Args, Parser, Subcommand};
use seqtk_rs::fq_check;

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// fastq QC
    Fqchk(FqchkArgs),
}

#[derive(Args)]
struct FqchkArgs {
    /// input fastq path
    in_fq: String,
    /// output tsv path
    out: String,

    #[arg(short, long)]
    /// quality value (default: 0)
    quality_value: Option<u8>,

    /// ascii value (default: 33)
    ascii_base: Option<u8>,
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Fqchk(fqchk) => {
            fq_check::fq_check(
                &fqchk.in_fq,
                &fqchk.out,
                fqchk.quality_value.unwrap_or(0),
                fqchk.ascii_base.unwrap_or(33),
            )?;
        }
    }
    Ok(())
}
