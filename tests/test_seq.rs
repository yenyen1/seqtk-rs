#[cfg(test)]
mod tests {
    use std::process::Command;
    use std::{fs, str};

    fn run_program_with_args(args: &[&str]) -> String {
        let output = Command::new("cargo")
            .arg("run")
            .args(args)
            .output()
            .expect("Failed to execute command");

        str::from_utf8(&output.stdout)
            .expect("Invalid UTF-8 output")
            .to_string()
    }

    #[test]
    fn test_seq_input_fastq_output_fastq() {
        // 01 - do nothing
        let args: Vec<&str> = "seq -I tests/data/chr.fastq".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/chr.fastq").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fq - 01");
        // 02 - mask
        let args: Vec<&str> =
            "seq -I tests/data/chr.fastq --drop-ambigous-seq --mask-char x --lowercases-to-char"
                .split_whitespace()
                .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fq_02.fastq").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fq - 02");
        // 03 - mask by qual
        let args: Vec<&str> = "seq -I tests/data/chr.fastq -1 --q-low 5 --mask-char x -M tests/data/chr.bed --reverse-complement".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fq_03.fastq").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fq - 03");
        // 04 - uppercases + mask by qual + line length
        let args: Vec<&str> =
            "seq -I tests/data/chr.fastq --q-high 20 --uppercases --line-len 20 --trim-header"
                .split_whitespace()
                .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fq_04.fastq").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fq - 04");
        // 05 - mask by bed + coplement
        let args: Vec<&str> = "seq -I tests/data/chr.fastq -M tests/data/chr.bed --mask-complement-region --fake-fastq-quality ^".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fq_05.fastq").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fq - 05");
        // 06 - ascii bases
        let args: Vec<&str> =
            "seq -I tests/data/long.fastq -l 1000 --ascii-bases 35 --output-qual-33"
                .split_whitespace()
                .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fq_06.fastq").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fq - 06");
    }
    #[test]
    fn test_seq_input_fastq_output_fasta() {
        // 01 - do nothing
        let args: Vec<&str> = "seq -I tests/data/chr.fastq --output-fasta"
            .split_whitespace()
            .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/chr.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fq - 01");
        // 02 -
        let args: Vec<&str> = "seq -I tests/data/chr.fastq --drop-ambigous-seq --mask-char x --lowercases-to-char --output-fasta".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fa_02.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fa - 02");
        // 03 -
        let args: Vec<&str> = "seq -I tests/data/chr.fastq -1 --q-low 5 --mask-char x -M tests/data/chr.bed --reverse-complement --output-fasta".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fa_03.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fa - 03");
        // 04 -
        let args: Vec<&str> = "seq -I tests/data/chr.fastq --q-high 20 --uppercases --line-len 20 --trim-header --output-fasta".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fa_04.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fa - 04");
        // 05 -
        let args: Vec<&str> = "seq -I tests/data/chr.fastq -M tests/data/chr.bed --mask-complement-region --fake-fastq-quality ^ --output-fasta".split_whitespace().collect();
        let output = run_program_with_args(&args);
        assert_eq!(output.as_bytes(), [], "[test] fq -> fa - 05");
        // 06 - ascii bases
        let args: Vec<&str> =
            "seq -I tests/data/long.fastq -l 1000 --ascii-bases 35 --q-low 5 --output-fasta"
                .split_whitespace()
                .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fa_06.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fq -> fa - 06");
    }
    #[test]
    fn test_seq_input_fasta_output_fasta() {
        // 01 - output fasta
        let args: Vec<&str> = "seq -A tests/data/chr.fasta --output-fasta"
            .split_whitespace()
            .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/chr.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fa - 01");
        // 02 -
        let args: Vec<&str> = "seq -A tests/data/chr.fasta --drop-ambigous-seq --mask-char x --lowercases-to-char --output-fasta".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fa_02.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fa - 02");
        // 03 -
        let args: Vec<&str> = "seq -A tests/data/chr.fasta -1 --mask-char x -M tests/data/chr.bed --reverse-complement --output-fasta".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fa2fa_03.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fa - 03");
        // 04 -
        let args: Vec<&str> =
            "seq -A tests/data/chr.fasta --uppercases --line-len 20 --trim-header --output-fasta"
                .split_whitespace()
                .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fa2fa_04.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fa - 04");
        // 05 -
        let args: Vec<&str> = "seq -A tests/data/chr.fasta -M tests/data/chr.bed --mask-complement-region --fake-fastq-quality ^ --output-fasta".split_whitespace().collect();
        let output = run_program_with_args(&args);
        assert_eq!(output.as_bytes(), [], "[test] fa -> fa - 05");
    }
    #[test]
    fn test_seq_input_fasta_output_fastq() {
        // 01 - output fastq
        let args: Vec<&str> = "seq -A tests/data/chr.fasta".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/chr.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fq - 01");
        // 02 -
        let args: Vec<&str> =
            "seq -A tests/data/chr.fasta --drop-ambigous-seq --mask-char x --lowercases-to-char"
                .split_whitespace()
                .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fa_02.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fq - 02");
        // 03 -
        let args: Vec<&str> = "seq -A tests/data/chr.fasta -1 --mask-char x -M tests/data/chr.bed --reverse-complement".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fa2fa_03.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fq - 03");
        // 04 -
        let args: Vec<&str> =
            "seq -A tests/data/chr.fasta --uppercases --line-len 20 --trim-header"
                .split_whitespace()
                .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fa2fa_04.fasta").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fq - 04");
        // 05 -
        let args: Vec<&str> = "seq -A tests/data/chr.fasta -M tests/data/chr.bed --mask-complement-region --fake-fastq-quality ^".split_whitespace().collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/seq/fq2fq_05.fastq").expect("");
        assert_eq!(output.as_bytes(), expect_content, "[test] fa -> fq - 05");
    }
}
