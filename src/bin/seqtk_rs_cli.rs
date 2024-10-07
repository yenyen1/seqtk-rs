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

    #[arg(short, long)]
    /// quality value (default: 0)
    quality_value: Option<u8>,
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Fqchk(fqchk) => {
            println!("'myapp add' was used, name is: {:?}", fqchk.in_fq);
            fq_check::fq_check(&fqchk.in_fq, fqchk.quality_value);
        }
    }
}
