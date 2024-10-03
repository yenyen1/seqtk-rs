use clap::{Args, Parser, Subcommand};
use seqtk_rs::abc;

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
    /// quality value
    quality_value: Option<u32>,
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Fqchk(fqchk) => {
            println!("'myapp add' was used, name is: {:?}", fqchk.in_fq);
            abc(&fqchk.in_fq, fqchk.quality_value);
        }
    }
}
