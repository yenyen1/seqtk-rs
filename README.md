# seqtk-rs
[![crate](https://img.shields.io/crates/v/seqtk-rs.svg)](https://crates.io/crates/seqtk-rs)

This is a sequence processing tool written in Rust for manipulating FASTA/FASTQ files. I built this tool out of my passion for Rust. Its functionality and subcommand names are similar to those in [`seqtk`](https://github.com/lh3/seqtk), but I’ve made some changes based on my own design logic. 

<!-- ⚠️ **Notice:** This project was previously paused but is now being actively resumed. The first release is expected by the end of June or early July 2025. -->

## Installation
```sh
# through cargo
cargo install seqtk-rs

# or through pip
pip install seqtk-rs

seqtk_rs -h
```

## Current Features
- [x] `seq`     Common transformation of FASTA/Q
                
- [x] `sample`  Random Sampling by given seed and fraction

- [x] `size`    Report the stats of sequence length 
  
    (**Output:** #seq, #bases, avg_size, min_size, med_size, max_size, N50)


- [x] `fqchk`   Report stats for sequence and quality by position
  
    (**Output:** POS, #bases, %A, %C, %G, %T, %N, avgQ, errQ, ...)
    - **avgQ:** Average quality score *`(Q₁ + Q₂ + ... + Qₙ) / N`*
    - **errQ:** Estimated error rate *`-10 * log₁₀((P₁ + P₂ + ... + Pₙ) / N)`*

    **Notice:** Some tools treat quality scores less than 3 (Q < 3) as 3 to avoid instability in downstream metrics. For example, Q = 0 yields an error probability P = 1.0, Q = 1 gives P ≈ 0.794, and Q = 2 gives P ≈ 0.630. These low Q-scores can heavily skew error rate calculations (e.g., errQ), which is why they are often floored to 3. However, this adjustment can lead to results that are inconsistent with the original definition. Therefore, this tool preserves the original quality scores as-is.
  
- [x] `comp`    Report the nucleotide composition of FASTA/Q 
    
    (**Output**: #A, #C, #G, #T, #2, #3, #4, #CG, #GC)

    - `CG` or `GC`: Number of CG/GC on the template strand

- [x] `qctrim`    Trims low-quality bases from a FASTQ data based on a quality threshold Q.
    

## TODO
- [ ] `trimAdapter` trim the adapter for FASTQ file 


## Acknowledgements
- [`seqtk`](https://github.com/lh3/seqtk)