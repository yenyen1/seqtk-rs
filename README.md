# seqtk-rs
This is a sequence processing tool written in Rust that replicates the core functionality of seqtk for manipulating FASTA/FASTQ files.

<!-- ⚠️ **Notice:** This project was previously paused but is now being actively resumed. The first release is expected by the end of June or early July 2025. -->

## Installation
```sh
cargo install seqtk-rs
seqtk_rs -h
```

## Current supported features
* `seq`     Common transformation of FASTA/Q
* `sample`  Random Sampling by given seed and fraction
* `size`    Report the stats of sequence length (Output: #seq, #bases, avg_size, min_size, med_size, max_size, N50)
* `fqchk`   Report stats for sequence and quality by position (Output: POS, #bases, %A, %C, %G, %T, %N, avgQ, errQ, ...)
* `comp`    Report the nucleotide composition of FASTA/Q (Output: #A, #C, #G, #T, #2, #3, #4, #CG, #GC)