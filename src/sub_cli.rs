use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Common transformation of FASTA/Q
    Seq(SeqArgs),
    /// fastq QC summary
    Fqchk(FqchkArgs),
    /// Random Sampling given seed and fraction
    Sample(SampleArgs),
}

#[derive(Args)]
pub struct FqchkArgs {
    /// input fastq path
    pub in_fq: String,
    #[arg(short, long)]
    /// output tsv path
    pub out: Option<String>,
    #[arg(short, long)]
    /// quality value [default: 0]
    pub quality_value: Option<u8>,
    #[arg(short, long)]
    /// ascii value [default: 33]
    pub ascii_base: Option<u8>,
}

#[derive(Args)]
pub struct SampleArgs {
    #[arg(short = 'I', long)]
    /// Input fastq path
    pub in_fq: Option<String>,
    #[arg(short = 'A', long)]
    /// Input fasta path
    pub in_fa: Option<String>,
    #[arg(short = 's', long)]
    /// Set the seed for the random number generator. This value ensures reproducibility of the sampling process. (This option takes effect only when used in conjunction with --sample-fraction / -f.) [default: 4]
    pub random_seed: Option<u64>,
    #[arg(short = 'f', long)]
    /// Specify the fraction of the total dataset to sample. The value is a FLOAT between 0 and 1. For example, a value of 0.1 will sample 10% of the data.
    pub sample_fraction: Option<f64>,
}

#[derive(Args)]
pub struct SeqArgs {
    #[arg(short = 'I', long)]
    /// Input fastq path
    pub in_fq: Option<String>,
    #[arg(short = 'A', long)]
    /// Input fasta path
    pub in_fa: Option<String>,

    #[arg(short = 'L', long)]
    /// Remove sequences shorter than MINI_SEQ_LENGTH. [default: 0]
    pub mini_seq_length: Option<usize>,
    #[arg(short = 'N', long)]
    /// drop sequences containing ambiguous bases 'N'
    pub drop_ambigous_seq: bool,
    #[arg(short = '1', long)]
    /// Output only the reads from odd-numbered (2n-1) records 
    pub output_odd: bool,
    #[arg(short = '2', long)]
    /// Output only the reads from even-numbered (2n) records 
    pub output_even: bool,

    #[arg(short = 'r', long)]
    /// reverse complement
    pub reverse_complement: bool,
    #[arg(short = 'R', long)]
    /// output both forward and reverse complement
    pub both_complement: bool,
    #[arg(long)]
    /// force output format to FASTA (discard quality)
    pub output_fasta: bool,
    #[arg(short= 'C', long)]
    /// drop comments at the header lines (only keep the first word before first space)
    pub trim_header: bool,
    #[arg(short = 'l', long)]
    /// Number of characters per line for sequences and their corresponding quality values [default: all on a single line]
    pub line_len: Option<usize>,

    #[arg(short = 'Q', long)]
    /// The quality scores are represented as chars with ASCII values equal to the score plus a base offset ASCII_BASES [default: 33]
    pub ascii_bases: Option<u8>,
    #[arg(long)]
    /// Output the quality score to an offset of 33 (Effective only when --ascii-bases is not 33)
    pub output_qual_33: bool,
    #[arg(long)]
    /// Mask bases with a quality score lower than Q_LOW [default: 0]
    pub q_low: Option<u8>,
    #[arg(long)]
    /// Mask bases with a quality score higher than Q_HIGH [default: 255]
    pub q_high: Option<u8>,
    #[arg(short = 'F', long)]
    /// Generate fake quality values using the specified CHAR.
    pub fake_fastq_quality: Option<char>,

    #[arg(short = 'U', long)]
    /// Converts all bases in the sequences to uppercase. When used in conjunction with other masking options
    /// (e.g., --q-low, --q-high, --mask-regions, --mask-char, etc.),
    /// the program first converts the sequences to uppercase and then applies the other masking operations.
    pub uppercases: bool,
    #[arg(short = 'x', long)]
    /// Convert all lowercases to --mark-char
    pub lowercases_to_char: bool,

    #[arg(long)]
    /// Mask bases by converting them to MASK_CHAR [default: convert to lowercase]
    pub mask_char: Option<char>,
    #[arg(short = 'M', long)]
    /// Mask bases that overlap with the regions specified in the BED (0-based) file [default: null]
    pub mask_regions: Option<String>,
    #[arg(long)]
    /// Mask bases that do NOT overlap with the region specified in the BED (effective with --mask-regions / -M)
    pub mask_complement_region: bool,
}
