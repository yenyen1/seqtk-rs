use clap::{Args, Parser, Subcommand};

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
    name: Option<String>,
}

fn main() {
    let cli = Cli::parse();

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Commands::Fqchk(name) => {
            println!("'myapp add' was used, name is: {:?}", name.name);
        }
    }
}
