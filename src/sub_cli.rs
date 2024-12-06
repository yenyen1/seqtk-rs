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
    /// fastq QC summary
    Fqchk(FqchkArgs),
    /// Seq
    Seq(SeqArgs),
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
pub struct SeqArgs {
    #[arg(short = 'i', long)]
    /// input fastq or fasta path
    pub in_fx: String,
    /// out path
    #[arg(short = 'o', long)]
    pub out: Option<String>,

    #[arg(long)]
    /// mask bases that quality lower than q_low [default: 0]
    q_low: Option<u8>,
    #[arg(long)]
    /// mask bases that quality higher than q_high [default: 255]
    q_high: Option<u8>,
    #[arg(long)]
    /// mask bases converted to CHAR [default: convert to lowercase]
    mask_char: Option<char>,
    #[arg(long)]
    /// number of residues per line; 0 for 2^32-1 [default: 0]
    n_residues: Option<u32>,
    #[arg(long)]
    /// quality shift: ASCII-INT gives base quality [default: 33]
    quality_shift: Option<u8>,
    #[arg(short = 's', long)]
    /// random seed (effective with --sample-fraction / -f) [default: 11]
    pub random_seed: Option<u64>,
    #[arg(short = 'f', long)]
    /// sample FLOAT fraction of sequences [default: 1.]
    pub sample_fraction: Option<f64>,
    #[arg(short = 'M', long)]
    /// mask regions in BED or name list FILE [default: null]
    pub mask_regions: Option<String>,
    #[arg(short = 'L', long)]
    /// drop sequences with length shorter than INT [default: 0]
    pub mini_seq_length: Option<usize>,
    #[arg(long)]
    /// fake FASTQ quality []
    fake_fastq_quality: bool,
    #[arg(long)]
    /// mask complement region (effective with --mask-regions / -M)
    pub mask_complement_region: bool,
    #[arg(long)]
    /// reverse complement
    reverse_complement: bool,
    #[arg(long)]
    /// output both forward and reverse complement
    both_complement: bool,
    #[arg(long)]
    /// force FASTA output (discard quality)
    output_fasta: bool,
    #[arg(long)]
    /// drop comments at the header lines
    trim_header: bool,
    #[arg(long)]
    /// drop sequences containing ambiguous bases
    pub drop_ambigous_seq: bool,
    #[arg(short = '1', long)]
    /// output the 2n-1 reads only
    pub output_odd_reads: bool,
    #[arg(short = '2', long)]
    /// output the 2n reads only
    pub output_even_reads: bool,
    #[arg(long)]
    /// shift quality by '(-Q) - 33'
    shift_quality_33: bool,
    #[arg(short = 'U', long)]
    /// convert all bases to uppercases
    uppercases: bool,
    #[arg(short = 'x', long)]
    /// convert all lowercases to -n
    lowercases: bool,
    #[arg(short = 'S', long)]
    /// strip of white spaces in sequences
    strip_whitespace: bool,
}
