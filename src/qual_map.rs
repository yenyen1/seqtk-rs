use crate::dna::DNA;
use crate::sequence_read::SequenceRead;
use std::collections::HashMap;

pub enum QualMap {
    /// [ Qual count map ]
    /// If qual threshold is zero, print out the number of qual count for each position.
    /// If qual threshold is not zero, print out the % of low/high qual for each position.
    DNACount(HashMap<usize, HashMap<DNA, usize>>), // {Pos: {DNA: count}}
    QualCount(HashMap<usize, HashMap<usize, usize>>), // {Qual: {Pos: count}}

    /// [ Low qual sum ]
    /// Calc avgQ and errQ for each position.
    QualSum(HashMap<char, HashMap<usize, f64>>), // {ScoreType: {Pos: sum of Qual}}
}

impl QualMap {
    /// [Initial fqchk map]
    /// Initial DNA count map {Pos: {DNA: count}}
    pub fn init_dna_count_map() -> QualMap {
        let map = QualMap::init_dna_count();
        QualMap::DNACount(map)
    }
    /// Initial Qual sum map {ScoreType: {Pos: sum of Qual}}
    pub fn init_qual_sum_map() -> QualMap {
        let map = HashMap::from([
            ('Q', HashMap::from([(0, 0.0)])), // save sum of Q score for each pos
            ('P', HashMap::from([(0, 0.0)])), // save sum of P_error for each pos
        ]);
        QualMap::QualSum(map)
    }
    /// Initial QualMap depend on the qual threshold.
    pub fn init_qual_count_map(qual_threshold: u8) -> QualMap {
        match qual_threshold {
            0 => {
                let map = QualMap::init_pos_qual_count();
                QualMap::QualCount(map)
            }
            _ => {
                let map = QualMap::init_dna_count();
                QualMap::DNACount(map)
            }
        }
    }
    fn init_pos_qual_count() -> HashMap<usize, HashMap<usize, usize>> {
        HashMap::from([(0, HashMap::from([(0, 0)]))])
    }
    fn init_dna_count() -> HashMap<usize, HashMap<DNA, usize>> {
        HashMap::from([(0, QualMap::init_inner_dna_count_map())]) // Use pos = 0 to restore total count
    }
    fn init_inner_dna_count_map() -> HashMap<DNA, usize> {
        let mut dna_hashmap = HashMap::new();
        for c in DNA::all_variants() {
            dna_hashmap.insert(c, 0);
        }
        dna_hashmap
    }

    /// [ Update data to map ]
    /// Update data with a sequence read
    pub fn update_map_with_read(qual_map: &mut QualMap, read: &SequenceRead, qual_threshold: u8) {
        match qual_map {
            QualMap::QualSum(map) => {
                QualMap::update_qual_sum(map, read);
            }
            QualMap::QualCount(map) => {
                QualMap::update_qual_count(map, read);
            }
            QualMap::DNACount(map) => {
                QualMap::update_dna_count(map, read, qual_threshold);
            }
        }
    }
    fn update_dna_count(
        dna_count_map: &mut HashMap<usize, HashMap<DNA, usize>>,
        read: &SequenceRead,
        qual_threshold: u8,
    ) {
        fn add_dna_count(
            pos: usize,
            dna_count_map: &mut HashMap<usize, HashMap<DNA, usize>>,
            char: &DNA,
        ) {
            let inner_map = dna_count_map
                .entry(pos)
                .or_insert_with(QualMap::init_inner_dna_count_map);
            *inner_map.entry(char.clone()).or_insert(0) += 1;
        }
        match qual_threshold {
            0 => {
                for (i, dna) in read.get_seq().iter().enumerate() {
                    add_dna_count(0, dna_count_map, dna);
                    add_dna_count(i + 1, dna_count_map, dna);
                }
            }
            _ => {
                // If qual threshold is not zero, we count dna with low qual
                for (i, qual) in read.get_q_score_vec().iter().enumerate() {
                    if *qual < qual_threshold {
                        let dna = read.get_seq_char(i);
                        add_dna_count(0, dna_count_map, &dna);
                        add_dna_count(i + 1, dna_count_map, &dna);
                    }
                }
            }
        }
    }
    fn update_qual_sum(qual_map: &mut HashMap<char, HashMap<usize, f64>>, read: &SequenceRead) {
        fn add(
            qual_map: &mut HashMap<char, HashMap<usize, f64>>,
            qual_vec: Vec<f64>,
            qual_type: char,
        ) {
            let inner_map = qual_map.get_mut(&qual_type).expect("map should exist");
            for (i, q) in qual_vec.iter().enumerate() {
                *inner_map.entry(0).or_insert(0.0) += q;
                *inner_map.entry(i + 1).or_insert(0.0) += q;
            }
        }
        let q_score_vec = read.get_q_score_vec_f64(); // Q score
        add(qual_map, q_score_vec, 'Q');
        let p_err_vec = read.get_p_err_vec(); // P_error
        add(qual_map, p_err_vec, 'P');
    }
    fn update_qual_count(
        qual_map: &mut HashMap<usize, HashMap<usize, usize>>,
        read: &SequenceRead,
    ) {
        // need to add
    }

    /// [Stats output]
    /// Get {%Low}\t{%High}\t{%A}\t{%C}\t{%G}\t{%T} String from low_qual_count_map
    pub fn get_low_qual_stats(low_qual_count_map: &QualMap, pos: usize, total: f64) -> String {
        let low_qual_total: usize = QualMap::sum_map_count(low_qual_count_map, pos);
        let p_low = 100.0 * low_qual_total as f64 / total;
        let out = format!(
            "{:.1}\t{:.1}\t{}",
            p_low,
            100.0 - p_low,
            QualMap::get_count_proportion_stats(low_qual_count_map, pos, low_qual_total as f64)
        );
        out
    }
    pub fn get_dna_count_stats(dna_count_map: &QualMap, pos: usize, total: f64) -> String {
        match pos {
            0 => {
                let out = format!(
                    "Total\t{}",
                    QualMap::get_count_proportion_stats(dna_count_map, pos, total)
                );
                out
            }
            _ => {
                let out = format!(
                    "{}\t{}",
                    pos,
                    QualMap::get_count_proportion_stats(dna_count_map, pos, total)
                );
                out
            }
        }
    }
    /// Get {avgQ}\t{Qerr} String from qual_sum_map
    pub fn get_qual_sum_stats(qual_sum_map: &QualMap, pos: usize, total: f64) -> String {
        fn get_sum_qual_from_map(qual_sum_map: &QualMap, mode: char, pos: usize) -> f64 {
            match qual_sum_map {
                QualMap::QualSum(map) => {
                    let inner_map = map.get(&mode).expect("map should exist");
                    *inner_map.get(&pos).expect("pos in Q map should exist")
                }
                _ => {
                    panic!("should not happen");
                } // should not happen
            }
        }
        let q_avg = get_sum_qual_from_map(qual_sum_map, 'Q', pos) / total;
        let p_avg = get_sum_qual_from_map(qual_sum_map, 'P', pos) / total;
        let q_err = SequenceRead::convert_p_err_to_q_score(p_avg);
        let out = format!("{}\t{}", q_avg, q_err);
        out
    }
    /// {%A}\t{%C}\t{%G}\t{%T}\t{%N} String from low_qual_count_map
    fn get_count_proportion_stats(qual_map: &QualMap, pos: usize, total: f64) -> String {
        match qual_map {
            QualMap::DNACount(map) => {
                fn get_prop(inner_map: Option<&HashMap<DNA, usize>>, dna: DNA, total: f64) -> f64 {
                    match inner_map {
                        None => 0.0,
                        Some(map) => {
                            let count = map.get(&dna).expect("should not happen");
                            *count as f64 / total
                        }
                    }
                }
                let inner_map = map.get(&pos);
                let mut out = String::new();
                for dna in DNA::all_variants() {
                    // GET A,C,G,T,N by order
                    let tmp_str = get_prop(inner_map, dna, total);
                    out.push_str(&format!("\t{:.1}", tmp_str));
                }
                out
            }
            QualMap::QualCount(map) => {
                let out = String::new();
                out
            }
            _ => panic!("should not happen"),
        }
    }

    /// [Other function]
    /// Calculate the sum of count from map
    pub fn sum_map_count(qual_map: &QualMap, pos: usize) -> usize {
        match qual_map {
            QualMap::DNACount(map) => {
                let inner_map = map.get(&pos);
                match inner_map {
                    None => 0,
                    Some(map) => map.values().sum(),
                }
            }
            QualMap::QualCount(map) => {
                let inner_map = map.get(&pos);
                match inner_map {
                    None => 0,
                    Some(map) => map.values().sum(),
                }
            }
            _ => panic!("should not happen"),
        }
    }
}
