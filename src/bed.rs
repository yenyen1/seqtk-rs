use rayon::prelude::*;
use std::collections::HashMap;

pub struct BedMap {
    map: HashMap<String, Vec<[usize; 2]>>,
}
impl BedMap {
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
        }
    }
    pub fn add(&mut self, key: String, region: [usize; 2]) {
        self.map.entry(key).or_default().push(region);
    }
    pub fn get(&self, key: &str) -> Option<&Vec<[usize; 2]>> {
        self.map.get(key)
    }
    /// Merge overlapping intervals in the map and return a new HashMap
    pub fn merge(&mut self) {
        let merged_map: HashMap<String, Vec<[usize; 2]>> = self
            .map
            .par_iter()
            .map(|(name, intervals)| {
                let mut intervals: Vec<[usize; 2]> = intervals.to_vec();
                if intervals.len() <= 1 {
                    return (name.clone(), intervals);
                }
                // Sort by start position
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
                (name.clone(), merged_intervals)
            })
            .collect();
        self.map = merged_map;
    }
}

pub fn is_overlapping(pos: usize, bed_pos: &[[usize; 2]]) -> (bool, usize) {
    // (is_overlap, end_idx)
    if bed_pos.is_empty() {
        (false, 0)
    } else {
        for (idx, &range) in bed_pos.iter().enumerate() {
            if range[1] < pos {
            } else if pos < range[1] && range[0] <= pos {
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
        for region in [[1, 5], [6, 10], [2, 7]] {
            bed_map.add("id_1".to_string(), region);
        }
        for region in [[13, 21], [31, 35], [19, 23], [6, 10], [30, 32], [41, 45]] {
            bed_map.add("id_2".to_string(), region);
        }
        bed_map.merge();
        let id_1_bed = bed_map.get("id_1").unwrap();
        let id_2_bed = bed_map.get("id_2").unwrap();
        assert_eq!(*id_1_bed, vec![[1, 10]]);
        assert_eq!(*id_2_bed, vec![[6, 10], [13, 23], [30, 35], [41, 45]]);
    }
    #[test]
    fn test_is_overlapping() {
        let bed_pos = vec![[6, 10], [13, 23], [30, 35], [41, 45]];
        let result = is_overlapping(11, &bed_pos[0..]);
        assert_eq!(result, (false, 1));
        let result = is_overlapping(31, &bed_pos[1..]);
        assert_eq!(result, (true, 1));
        let result = is_overlapping(32, &bed_pos[2..]);
        assert_eq!(result, (true, 0));
    }
}
