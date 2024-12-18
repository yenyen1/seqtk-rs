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
    #[arg(short = 'I', long)]
    /// input fastq path
    pub in_fq: Option<String>,
    #[arg(short = 'A', long)]
    /// input fasta path
    pub in_fa: Option<String>,
    #[arg(long)]
    /// ascii bases [default: 33]
    pub ascii_bases: Option<u8>,
    #[arg(long)]
    /// mask bases that quality lower than q_low [default: 0]
    pub q_low: Option<u8>,
    #[arg(long)]
    /// mask bases that quality higher than q_high [default: 255]
    pub q_high: Option<u8>,
    #[arg(long)]
    /// mask bases converted to CHAR [default: convert to lowercase]
    pub mask_char: Option<char>,
    #[arg(long)]
    /// number of characters per line for seqence and quality [default: all in one line]
    pub line_len: Option<usize>,
    #[arg(short = 's', long)]
    /// random seed (effective with --sample-fraction / -f) [default: 4]
    pub random_seed: Option<u64>,
    #[arg(short = 'f', long)]
    /// sample FLOAT fraction of sequences [default: 1.]
    pub sample_fraction: Option<f64>,
    #[arg(short = 'M', long)]
    /// mask regions in BED or name list FILE [default: null]
    pub mask_regions: Option<String>,
    #[arg(short = 'l', long)]
    /// drop sequences with length shorter than INT [default: 0]
    pub mini_seq_length: Option<usize>,
    #[arg(long)]
    /// output fake FASTQ quality char
    pub fake_fastq_quality: Option<char>,
    #[arg(long)]
    /// mask complement region (effective with --mask-regions / -M)
    pub mask_complement_region: bool,
    #[arg(short = 'r', long)]
    /// reverse complement
    pub reverse_complement: bool,
    #[arg(short = 'R', long)]
    /// output both forward and reverse complement
    pub both_complement: bool,
    #[arg(long)]
    /// force FASTA output (discard quality)
    pub output_fasta: bool,
    #[arg(long)]
    /// drop comments at the header lines
    pub trim_header: bool,
    #[arg(short = 'N', long)]
    /// drop sequences containing ambiguous bases
    pub drop_ambigous_seq: bool,
    #[arg(short = '1', long)]
    /// output the 2n-1 reads only
    pub output_odd_reads: bool,
    #[arg(short = '2', long)]
    /// output the 2n reads only
    pub output_even_reads: bool,
    #[arg(long)]
    /// output by quality that shift with value 'Q'.
    pub output_shift_qual: Option<u8>,
    #[arg(short = 'U', long)]
    /// convert all bases to uppercases
    pub uppercases: bool,
    #[arg(short = 'x', long)]
    /// convert all lowercases to -n
    pub lowercases_to_char: bool,
    #[arg(short = 'S', long)]
    /// strip of white spaces in sequences
    pub strip_whitespace: bool,
}
