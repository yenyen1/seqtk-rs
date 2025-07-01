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
    fn test_fqchk() {
        // 01 - one record
        let args: Vec<&str> = "fqchk tests/data/test_cp.fastq"
            .split_whitespace()
            .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/fqchk/result_for_test_cp.txt").expect("");
        assert_eq!(output.as_bytes(), expect_content, "Err1");

        
        // 02 - two recrod
        let args: Vec<&str> = "fqchk  tests/data/fqchk/test_fqchk.fastq -q 10"
            .split_whitespace()
            .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/fqchk/result_for_test_fqchk.txt").expect("");
        assert_eq!(output.as_bytes(), expect_content, "Err2");

        // 03 - test Q < 3
        let args: Vec<&str> = "fqchk -q 10 tests/data/fqchk/test_fqchk2.fastq"
            .split_whitespace()
            .collect();
        let output = run_program_with_args(&args);
        let expect_content = fs::read("tests/data/fqchk/result_for_test_fqchk2.txt").expect("");
        assert_eq!(output.as_bytes(), expect_content, "Err3");

    }
}
