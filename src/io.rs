use crate::bed::BedMap;

use std::fs::File;
use std::io::{self, BufRead, BufReader};

pub fn get_bed_map(file_path: &str) -> Result<BedMap, io::Error> {
    // Read the file and collect unique intervals for each name
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut bed_map: BedMap = BedMap::new();
    for line in reader.lines() {
        match line {
            Ok(line_content) => {
                let columns: Vec<&str> = line_content.split('\t').collect();
                if columns.len() > 2 {
                    let name = columns[0].to_string();
                    if let (Ok(start), Ok(end)) =
                        (columns[1].parse::<usize>(), columns[2].parse::<usize>())
                    {
                        bed_map.add(name, [start, end]);
                    } else {
                        eprintln!("Error parsing start or end for line: {}", line_content);
                    }
                }
            }
            Err(e) => eprintln!("Error reading line: {}", e),
        }
    }
    bed_map.merge();
    Ok(bed_map)
}
