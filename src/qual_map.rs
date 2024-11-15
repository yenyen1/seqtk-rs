use crate::dna::DNA;
use crate::sequence_read::SequenceRead;
use std::collections::HashMap;

pub enum QualMap {
    /// [ Qual count map ]
    /// If qual threshold is zero, print out the number of qual count for each position.
    /// If qual threshold is not zero, print out the % of low/high qual for each position.
    DNACount(HashMap<usize, HashMap<DNA, usize>>), // {Pos: {DNA: count}}
    QualCount(HashMap<usize, HashMap<usize, usize>>), // {Qual: {Pos: count}}
    LowQualCount(HashMap<usize,usize>), // {Pos: low_qual_count}

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
                let map = QualMap::init_low_qual_count();
                QualMap::LowQualCount(map)
            }
        }
    }
    fn init_pos_qual_count() -> HashMap<usize, HashMap<usize, usize>> {
        HashMap::from([(0, HashMap::from([(0, 0)]))])
    }
    fn init_low_qual_count() -> HashMap<usize,usize> {
        HashMap::from([(0,0)])
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
                QualMap::update_dna_count(map, read);
            }
            QualMap::LowQualCount(map) => {
                QualMap::update_low_qual_count(map, read, qual_threshold);
            }
        }
    }
    fn update_dna_count(
        dna_count_map: &mut HashMap<usize, HashMap<DNA, usize>>,
        read: &SequenceRead,
    ) {
        fn add_dna_count(
            pos: usize,
            dna_count_map: &mut HashMap<usize, HashMap<DNA, usize>>,
            dna: &DNA,
        ) {
            let inner_map = dna_count_map
                .entry(pos)
                .or_insert_with(QualMap::init_inner_dna_count_map);
            *inner_map.entry(dna.clone()).or_insert(0) += 1;
        }
        for (i, dna) in read.get_seq().as_slice().iter().enumerate() {
            add_dna_count(0, dna_count_map, dna);
            add_dna_count(i + 1, dna_count_map, dna);
        }
    }
    fn update_qual_sum(qual_map: &mut HashMap<char, HashMap<usize, f64>>, read: &SequenceRead) {
        fn add(
            qual_map: &mut HashMap<char, HashMap<usize, f64>>,
            qual_vec: &Vec<f64>,
            qual_type: char,
        ) {
            let inner_map = qual_map.get_mut(&qual_type).expect("0 map should exist");
            for (i, q) in qual_vec.iter().enumerate() {
                *inner_map.entry(0).or_insert(0.0) += q;
                *inner_map.entry(i + 1).or_insert(0.0) += q;
            }
        }
        let q_score_vec = read.get_q_score_vec_f64(); // Q score
        add(qual_map, &q_score_vec, 'Q');
        let p_err_vec = read.get_p_err_vec(); // P_error
        add(qual_map, &p_err_vec, 'P');
    }
    fn update_low_qual_count(qual_map: &mut HashMap<usize, usize>, read: &SequenceRead, qual_threshold: u8) {
        for (i, q) in read.get_q_score_vec().iter().enumerate() {
            if *q < qual_threshold {
                *qual_map.entry(0).or_insert(0) += 1;
                *qual_map.entry(i+1).or_insert(0) += 1;
            }
        }   
    }
    fn update_qual_count(
        qual_map: &mut HashMap<usize, HashMap<usize, usize>>,
        read: &SequenceRead,
    ) {
        // need to add
    }

    /// [Stats output]
    pub fn get_qual_count_stats(qual_map: &QualMap, pos: usize, total: f64) -> String {
        match qual_map {
            QualMap::LowQualCount(map) => {
                QualMap::get_low_qual_stats(map, pos, total)
            }
            QualMap::QualCount(_) => {
                String::new()
            }
            _ => panic!("1 should not happen")
        }
    }
    /// Get \t{%Low}\t{%High} String from low_qual_count_map
    fn get_low_qual_stats(low_qual_count_map: &HashMap<usize, usize>, pos: usize, total: f64) -> String {
        let low_qual_total = low_qual_count_map.get(&pos);
        let p_low = match low_qual_total {
            Some(&v) => 100.0 * v as f64 / total,
            None => 0.0,
        };
        let out = format!(
            "\t{:.1}\t{:.1}",
            p_low,
            100.0 - p_low,
        );
        out
    }
    pub fn get_dna_count_stats(dna_count_map: &QualMap, pos: usize, total: f64) -> String {
        let count_prop = QualMap::get_count_proportion_stats(dna_count_map, pos, total);
        match dna_count_map {
            QualMap::DNACount(_) => {
                match pos {
                    0 => format!("Total\t{:.0}{}",total,count_prop),
                    _ => format!("{}\t{:.0}{}",pos,total,count_prop),
                }
            }
            QualMap::QualCount(_) => {
                format!("{}",count_prop)
            }
            _ => panic!("2 should not happen")
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
                    panic!("3 should not happen");
                } // should not happen
            }
        }
        let q_avg = get_sum_qual_from_map(qual_sum_map, 'Q', pos) / total;
        let p_avg = get_sum_qual_from_map(qual_sum_map, 'P', pos) / total;
        let q_err = SequenceRead::convert_p_err_to_q_score(p_avg);
        format!("\t{:.1}\t{:.1}", q_avg, q_err)
    }
    /// {%A}\t{%C}\t{%G}\t{%T}\t{%N} String from dna_count_map
    /// 
    fn get_count_proportion_stats(qual_map: &QualMap, pos: usize, total: f64) -> String {
        match qual_map {
            QualMap::DNACount(map) => {
                fn get_prop(inner_map: Option<&HashMap<DNA, usize>>, dna: DNA, total: f64) -> f64 {
                    match inner_map {
                        None => 0.0,
                        Some(map) => {
                            let count = map.get(&dna).expect("4 should not happen");
                            100.0 * (*count as f64) / total
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
            _ => panic!("5 should not happen"),
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
            _ => panic!("6 should not happen"),
        }
    }
    pub fn get_pos_length(qual_map: &QualMap) -> usize {
        match qual_map {
            QualMap::DNACount(map) => {
                let length = map.keys().len();
                length
            }
            _ => panic!("7 should not happen"),
        }
    }
}
