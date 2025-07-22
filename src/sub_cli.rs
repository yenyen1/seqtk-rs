use clap::{ArgGroup, Args, Parser, Subcommand};
use colored::*;

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Common transformation of FASTA/Q")]
    Seq(SeqArgs),

    #[command(about = "Random Sampling by given seed and fraction")]
    Sample(SampleArgs),

    #[command(
        about = "Report the stats of sequence length (Output: #seq, #bases, avg_size, min_size, med_size, max_size, N50)"
    )]
    Size(SizeArgs),

    #[command(
        about = "Report stats for sequence and quality by position (Output: POS, #bases, %A, %C, %G, %T, %N, avgQ, errQ, ...)",
        long_about = "\x1b[1mFqchk\n\x1b[0m\
                    Parses FASTQ data with quality threshold and computes per-position statistics.\n\n\
                    \x1b[1;4mOutput columns:\x1b[0m\n\
                    (1) POS: Position in the read\n\
                    (2) #bases: Number of bases at this position\n\
                    (3) %A, %C, %G, %T, %N: Percentage of each nucleotide\n\
                    (4) avgQ: Average quality score (Q₁ + Q₂ + ... + Qₙ) / N\n\
                    (5) errQ: Estimated error rate -10 * log₁₀((P₁ + P₂ + ... + Pₙ) / N)\n\
                    (6-7) %low, %high: Percentage of the nucleotide that the quality scores below or above the threshold, respectively. (When q_threshold > 0)\n\
                    (6-) %QX: Percentage of the nucleotide that the quality scores is X (When q_threshold = 0)\n\n\
                    \x1b[1;4mNote:\x1b[0m\n\
                    Some tools treat quality scores less than 3 (Q < 3) as 3 to avoid instability in downstream metrics. \
                    For example, Q = 0 yields an error probability P = 1.0, Q = 1 gives P ≈ 0.794, and Q = 2 gives P ≈ 0.630. \
                    These low Q-scores can heavily skew error rate calculations (e.g., errQ), which is why they are often floored to 3. \
                    However, this adjustment can lead to results that are inconsistent with the original definition. \
                    Therefore, this tool preserves the original quality scores as-is."
    )]
    Fqchk(FqchkArgs),

    #[command(
        about = "Report the nucleotide composition of FASTA/Q (Output: #A, #C, #G, #T, #2, #3, #4, #CG, #GC)",
        long_about = "\x1b[1mComp\n\x1b[0m\
                    Report the nucleotide composition of FASTA/Q w/o masked sequences with a BED file.\n\n\
                    \x1b[1;4mOutput columns:\x1b[0m\n\
                    (1) #seq: Number of reads\n\
                    (2) #bases: Number of bases\n\
                    (3-6) #A, #C, #G, #T: Number of each nucleotide\n\
                    (7) #2: Number of R, Y, S, W, K, M\n\
                    (8) #3: Number of B, D, H, V\n\
                    (9) #4: Number of N\n\
                    (10) #CG: Number of CG on the template strand\n\
                    (11) #GC: Number of GC on the template strand"
    )]
    Comp(CompArgs),

    #[command(
        about = "Trims low-quality bases from a FASTQ data based on a quality threshold Q.",
        long_about = "\x1b[1mQCTrim\n\x1b[0m\
                    Trims low-quality bases from a FASTQ data based on a quality threshold Q.\n\n\
                    \x1b[1;4mThe algorithm:\n\x1b[0m\
                    (1) Scan the read from 5’ to 3’ to find the first base with quality >= Q. This is the starting position of the trimmed read.\n\
                    (2) Compute a running score from 5’ to 3’ \n\tscore(i) = score(i-1) + quality(i) - Q \n    where the initial score is score(start) = quality(start) - Q.\n\
                    (3) Find out the maximun score. This is the end of the trimmed read.\n\
                    (4) If all base quality are less than Q, the read is discarded.\n\
                    (5) If the trimmed read is shorter than min_len, we first extend it from the 3’ end. If still too short, extend form the 5’ end until min_len is reached or no more bases are available.\
                    \n\n\x1b[1;4mNotes:\n\x1b[0m\
                    Quality trimming is no longer necessary in most modern sequencing pipelines. Its usefulness depends on the sequencing technology and the goals of your downstream analysis."
    )]
    Qctrim(QCTrimArgs),
}

#[derive(Args)]
pub struct FqchkArgs {
    /// FASTQ path
    pub in_fq: String,
    #[arg(short, long)]
    /// Quality value [default: 0]
    pub quality_value: Option<u8>,
    #[arg(short, long)]
    /// Ascii value [default: 33]
    pub ascii_base: Option<u8>,
}

#[derive(Args)]
pub struct QCTrimArgs {
    /// FASTQ path
    pub in_fq: String,
    #[arg(short, long)]
    /// Quality threshold [default: 13]
    pub q_thershold: Option<u8>,
    #[arg(short, long)]
    /// Minimum allowed length for a read after trimming [default: 30]
    pub min_length: Option<usize>,
    #[arg(short, long)]
    /// Ascii value [default: 33]
    pub ascii_base: Option<u8>,
}

#[derive(Args)]
#[command(group(
    ArgGroup::new("exclusive_group")
        .args(["in_fq", "in_fa"])
        .required(true)
        .multiple(false)
))]
pub struct SizeArgs {
    #[arg(short = 'I', long)]
    /// FASTQ path
    pub in_fq: Option<String>,
    #[arg(short = 'A', long)]
    /// FASTA path
    pub in_fa: Option<String>,
}

#[derive(Args)]
#[command(group(
    ArgGroup::new("exclusive_group")
        .args(["in_fq", "in_fa"])
        .required(true)
        .multiple(false)
))]
pub struct CompArgs {
    #[arg(short = 'I', long)]
    /// FASTQ path
    pub in_fq: Option<String>,
    #[arg(short = 'A', long)]
    /// FASTA path
    pub in_fa: Option<String>,
    #[arg(short = 'u', long)]
    /// Only report unmasked bases [default: false]
    pub exclude_masked: bool,
    #[arg(short = 'r', long)]
    /// Report bases that overlap with the regions specified in the BED (0-based) file [default: null]
    pub in_bed: Option<String>,
}

#[derive(Args)]
#[command(group(
    ArgGroup::new("exclusive_group")
        .args(["in_fq", "in_fa"])
        .required(true)
        .multiple(false)
))]
pub struct SampleArgs {
    #[arg(short = 'I', long)]
    /// FASTQ path
    pub in_fq: Option<String>,
    #[arg(short = 'A', long)]
    /// FASTA path
    pub in_fa: Option<String>,
    #[arg(short = 's', long)]
    /// Set the seed for the random number generator. This value ensures reproducibility of the sampling process. (This option takes effect only when used in conjunction with --sample-fraction / -f.) [default: 4]
    pub random_seed: Option<usize>,
    #[arg(short = 'f', long, value_parser = validate_ratio)]
    /// Specify the fraction of the total dataset to sample. The value is a FLOAT between 0 and 1. For example, a value of 0.1 will sample 10% of the data.
    pub sample_fraction: Option<f64>,
}

#[derive(Args)]
#[command(group(
    ArgGroup::new("exclusive_group")
        .args(["in_fq", "in_fa"])
        .required(true)
        .multiple(false)
))]
pub struct SeqArgs {
    #[arg(short = 'I', long)]
    /// FASTQ path
    pub in_fq: Option<String>,
    #[arg(short = 'A', long)]
    /// FASTA path
    pub in_fa: Option<String>,

    #[arg(short = 'L', long)]
    /// Remove sequences shorter than MINI_SEQ_LENGTH [default: 0]
    pub mini_seq_length: Option<usize>,
    #[arg(short = 'N', long)]
    /// Drop sequences containing ambiguous bases 'N' [default: false]
    pub drop_ambigous_seq: bool,
    #[arg(short = '1', long)]
    /// Output only the reads from odd-numbered (2n-1) records
    pub output_odd: bool,
    #[arg(short = '2', long)]
    /// Output only the reads from even-numbered (2n) records
    pub output_even: bool,

    #[arg(short = 'r', long)]
    /// Reverse complement [default: false]
    pub reverse_complement: bool,
    #[arg(short = 'R', long)]
    /// Output both forward and reverse complement [default: false]
    pub both_complement: bool,
    #[arg(long)]
    /// Force output format to FASTA (discard quality)
    pub output_fasta: bool,
    #[arg(short = 'C', long)]
    /// Drop comments at the header lines (only keep the first word before first space)
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
/// Validate seq arguments.
pub fn valiation_seq_args(args: &SeqArgs) -> Result<(), std::io::Error> {
    let mut errors = Vec::new();
    if args.output_even && args.output_odd {
        errors.push("--output-even-reads and --output-odd-reads can not be used together.");
    }
    if args.mask_complement_region && args.mask_regions.is_none() {
        errors.push("--mask-complment-region requires --mask-regions.");
    }
    if args.lowercases_to_char && args.mask_char.is_none() {
        errors.push("--lowercases-to-char requires --mask-char.");
    }
    if args.output_fasta && (args.output_qual_33 || args.fake_fastq_quality.is_some()) {
        errors
            .push("--output-fasta can not be used with --output-qual-33 or --fake-fastq-quality.");
    }
    if args.output_qual_33 && args.fake_fastq_quality.is_some() {
        errors.push("--output-qual-33 and --fake-fastq-quality can not be used together.");
    }
    if args.reverse_complement && args.both_complement {
        errors.push("--reverse-complement and --both-complement can not be used together.");
    }
    if !errors.is_empty() {
        for error in errors {
            eprintln!("{} {}", "error:".red().bold(), error);
        }
        std::process::exit(1);
    }
    Ok(())
}

fn validate_ratio(s: &str) -> Result<f64, String> {
    let val: f64 = s
        .parse()
        .map_err(|_| "Must be a valid floating-point number".to_string())?;
    if (0.0..=1.0).contains(&val) {
        Ok(val)
    } else {
        Err("Value must be between 0.0 and 1.0 (inclusive)".to_string())
    }
}
