use bio::io::{fasta, fastq};
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Stdout};
use std::path::Path;

pub fn append_bufwriter(out_path: &str) -> io::Result<BufWriter<File>> {
    let file = OpenOptions::new()
        .create(true) // Create the file if it doesn't exist
        .append(true) // Open the file in append mode (no truncation)
        .open(out_path)?; // Open the file at the provided path
    let writer = BufWriter::new(file);
    Ok(writer)
}
pub fn new_fq_iterator(file_path: &str) -> io::Result<fastq::Reader<BufReader<Box<dyn BufRead>>>> {
    let file_extension = Path::new(file_path).extension().and_then(|s| s.to_str());
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = match file_extension {
        Some("gz") => {
            let gz_decoder = GzDecoder::new(file);
            Box::new(BufReader::new(gz_decoder))
        }
        _ => Box::new(BufReader::new(file)),
    };
    Ok(fastq::Reader::new(reader))
}
pub fn new_fa_iterator(file_path: &str) -> io::Result<fasta::Reader<BufReader<Box<dyn BufRead>>>> {
    let file_extension = Path::new(file_path).extension().and_then(|s| s.to_str());
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = match file_extension {
        Some("gz") => {
            let gz_decoder = GzDecoder::new(file);
            Box::new(BufReader::new(gz_decoder))
        }
        _ => Box::new(BufReader::new(file)),
    };
    Ok(fasta::Reader::new(reader))
}
pub trait RecordType {
    fn seq(&self) -> &[u8];
    fn id(&self) -> &str;
    fn desc(&self) -> Option<&str>;
    fn qual(&self) -> &[u8];
}
impl RecordType for fasta::Record {
    fn seq(&self) -> &[u8] {
        self.seq()
    }
    fn id(&self) -> &str {
        self.id()
    }
    fn desc(&self) -> Option<&str> {
        self.desc()
    }
    fn qual(&self) -> &[u8] {
        panic!("[RecordType::Fasta] no qual function.")
    }
}
impl RecordType for fastq::Record {
    fn seq(&self) -> &[u8] {
        self.seq()
    }
    fn id(&self) -> &str {
        self.id()
    }
    fn desc(&self) -> Option<&str> {
        self.desc()
    }
    fn qual(&self) -> &[u8] {
        self.qual()
    }
}
pub enum FxWriter {
    Fasta(fasta::Writer<Stdout>),
    Fastq(fastq::Writer<Stdout>),
}
impl FxWriter {
    pub fn new(is_fasta: bool) -> FxWriter {
        let std_out = io::stdout();
        if is_fasta {
            FxWriter::Fasta(fasta::Writer::new(std_out))
        } else {
            FxWriter::Fastq(fastq::Writer::new(std_out))
        }
    }
    pub fn write(
        &mut self,
        id: &str,
        seq: &[u8],
        desc: Option<&str>,
        qual: &[u8],
    ) -> io::Result<()> {
        match self {
            FxWriter::Fasta(ref mut out) => {
                out.write(id, desc, seq)?;
            }
            FxWriter::Fastq(ref mut out) => {
                out.write(id, desc, seq, qual)?;
            }
        }
        Ok(())
    }
}
pub fn get_bed_map(file_path: &str) -> Result<HashMap<String, Vec<[usize; 2]>>, io::Error> {
    // Read the file and collect unique intervals for each name
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut bed_map: HashMap<String, HashSet<[usize; 2]>> = HashMap::new();
    for line in reader.lines() {
        match line {
            Ok(line_content) => {
                let columns: Vec<&str> = line_content.split('\t').collect();
                if columns.len() > 2 {
                    let name = columns[0].to_string();
                    if let (Ok(start), Ok(end)) =
                        (columns[1].parse::<usize>(), columns[2].parse::<usize>())
                    {
                        bed_map
                            .entry(name)
                            .or_insert_with(HashSet::new)
                            .insert([start, end]);
                    } else {
                        eprintln!("Error parsing start or end for line: {}", line_content);
                    }
                }
            }
            Err(e) => eprintln!("Error reading line: {}", e),
        }
    }
    // Merge overlapping intervals and store the results in a Vec
    let result_map = merge_bed(&bed_map);
    Ok(result_map)
}

fn merge_bed(bed_map: &HashMap<String, HashSet<[usize; 2]>>) -> HashMap<String, Vec<[usize; 2]>> {
    let result_map: HashMap<String, Vec<[usize; 2]>> = bed_map
        .par_iter()
        .map(|(name, intervals)| {
            let mut intervals: Vec<[usize; 2]> = intervals.iter().cloned().collect();
            if intervals.len() == 1 {
                return (name.to_string(), intervals);
            } else {
                intervals.sort_unstable_by(|a, b| a[0].cmp(&b[0]));
                let mut merged_intervals = Vec::new();
                for i in 1..intervals.len() {
                    if intervals[i - 1][1] < intervals[i][0] {
                        merged_intervals.push(intervals[i - 1]);
                    } else {
                        intervals[i][0] = intervals[i - 1][0];
                    }
                }
                merged_intervals.push(intervals[intervals.len() - 1]);
                return (name.to_string(), merged_intervals);
            }
        })
        .collect();
    result_map
}

#[cfg(test)]
mod tests {
    use crate::utils::merge_bed;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn test_merge_bed() {
        let mut bed_map = HashMap::new();
        bed_map.insert("id_1".to_string(), HashSet::from([[1, 5], [6, 10], [2, 7]]));
        bed_map.insert(
            "id_2".to_string(),
            HashSet::from([[13, 21], [31, 35], [19, 23], [6, 10], [30, 32], [41, 45]]),
        );
        println!("{:?}", bed_map);
        let result_map = merge_bed(&bed_map);
        println!("{:?}", result_map);
        let id_1_bed = result_map.get("id_1").unwrap();
        let id_2_bed = result_map.get("id_2").unwrap();
        assert_eq!(*id_1_bed, vec![[1, 10]]);
        assert_eq!(*id_2_bed, vec![[6, 10], [13, 23], [30, 35], [41, 45]]);
    }
}
