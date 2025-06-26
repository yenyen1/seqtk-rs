use crate::io_utils::buffer_reader_maybe_gz;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::BufRead;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct BedPos(pub usize, pub usize); // (start, end)
pub struct BedMap {
    map: HashMap<String, Vec<BedPos>>, // {rname: vec[Bedpos..]}
}
impl BedMap {
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
        }
    }
    pub fn add(&mut self, key: String, region: BedPos) {
        self.map.entry(key).or_default().push(region);
    }
    pub fn get(&self, key: &str) -> Option<&Vec<BedPos>> {
        self.map.get(key)
    }
    pub fn from(path: &str) -> Result<BedMap, std::io::Error> {
        let reader = buffer_reader_maybe_gz(path)?;
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
                            bed_map.add(name, BedPos(start, end));
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
    /// Merge overlapping intervals in the map and return a new HashMap
    fn merge(&mut self) {
        let merged_map: HashMap<String, Vec<BedPos>> = self
            .map
            .par_iter()
            .map(|(name, intervals)| {
                let mut intervals: Vec<BedPos> = intervals.to_vec();
                if intervals.len() <= 1 {
                    return (name.clone(), intervals);
                }
                // Sort by start position
                intervals.sort_unstable_by(|a, b| a.0.cmp(&b.0));

                let mut merged_intervals = Vec::new();
                for i in 1..intervals.len() {
                    if intervals[i - 1].1 < intervals[i].0 {
                        merged_intervals.push(intervals[i - 1]);
                    } else {
                        intervals[i].0 = intervals[i - 1].0;
                    }
                }
                merged_intervals.push(intervals[intervals.len() - 1]);
                (name.clone(), merged_intervals)
            })
            .collect();
        self.map = merged_map;
    }
    
}

pub fn is_overlapping(pos: usize, bed_pos: &[BedPos]) -> (bool, usize) {
    // (is_overlap, end_idx)
    if bed_pos.is_empty() {
        (false, 0)
    } else {
        for (idx, &range) in bed_pos.iter().enumerate() {
            if range.1 < pos {
            } else if pos < range.1 && range.0 <= pos {
                return (true, idx); // If there's an overlap, return the index of the interval
            } else {
                return (false, idx);
            }
        }
        (false, bed_pos.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_bed() {
        let mut bed_map = BedMap::new();
        for [start, end] in [[1, 5], [6, 10], [2, 7]] {
            bed_map.add("id_1".to_string(), BedPos(start, end));
        }
        for [start, end] in [[13, 21], [31, 35], [19, 23], [6, 10], [30, 32], [41, 45]] {
            bed_map.add("id_2".to_string(), BedPos(start, end));
        }
        bed_map.merge();
        let id_1_bed = bed_map.get("id_1").unwrap();
        let id_2_bed = bed_map.get("id_2").unwrap();
        assert_eq!(*id_1_bed, vec![BedPos(1, 10)]);
        let result = vec![
            BedPos(6, 10),
            BedPos(13, 23),
            BedPos(30, 35),
            BedPos(41, 45),
        ];
        assert_eq!(*id_2_bed, result);
    }
    #[test]
    fn test_is_overlapping() {
        let bed_pos = vec![
            BedPos(6, 10),
            BedPos(13, 23),
            BedPos(30, 35),
            BedPos(41, 45),
        ];
        let result = is_overlapping(11, &bed_pos[0..]);
        assert_eq!(result, (false, 1));
        let result = is_overlapping(31, &bed_pos[1..]);
        assert_eq!(result, (true, 1));
        let result = is_overlapping(32, &bed_pos[2..]);
        assert_eq!(result, (true, 0));
    }
}
