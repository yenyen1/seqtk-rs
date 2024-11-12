use crate::dna::DNA;
use crate::sequence_read::SequenceRead;
use std::collections::HashMap;

pub enum QualMap {
    /// [ Qual count map ]
    /// If qual threshold is zero, print out the number of qual count for each position.
    /// If qual threshold is not zero, print out the % of low/high qual for each position.
    QualThresNonZero(HashMap<usize, HashMap<DNA, usize>>), // {Pos: {DNA: count}}
    QualThresZero(HashMap<usize, HashMap<usize, usize>>), // {Qual: {Pos: count}}

    /// [ Low qual sum ]
    /// Calc avgQ and errQ for each position.
    QualSum(HashMap<char, HashMap<usize, f64>>), // {ScoreType: {Pos: sum of Qual}}
}

impl QualMap {
    /// Initial QualMap depend on the qual threshold.
    pub fn init_qual_count_hashmap(qual_threshold: u8) -> QualMap {
        match qual_threshold {
            0 => {
                let map = QualMap::init_pos_qual_count();
                QualMap::QualThresZero(map)
            }
            _ => {
                let map = QualMap::init_low_qual_count();
                QualMap::QualThresNonZero(map)
            }
        }
    }
    fn init_pos_qual_count() -> HashMap<usize, HashMap<usize, usize>> {
        HashMap::from([(0, HashMap::from([(0, 0)]))])
    }
    fn init_low_qual_count() -> HashMap<usize, HashMap<DNA, usize>> {
        HashMap::from([(0, DNA::create_a_dna_count_map())]) // Use pos = 0 to restore total count
    }

    /// Initial QualSumMap
    pub fn init_qual_sum_map() -> QualMap {
        let map = HashMap::from([
            ('Q', HashMap::from([(0, 0.0)])), // save sum of Q score for each pos
            ('P', HashMap::from([(0, 0.0)])),
        ]); // save sum of P_error for each pos
        QualMap::QualSum(map)
    }

    /// get

    /// Update data with sequence read
    pub fn update_map_by_seqread(qual_map: &mut QualMap, read: &SequenceRead, qual_threshold: u8) {
        match qual_map {
            QualMap::QualSum(map) => {
                QualMap::update_qual_sum(map, read);
            }
            QualMap::QualThresZero(map) => {}
            QualMap::QualThresNonZero(map) => {
                QualMap::update_low_qual_count(map, read, qual_threshold);
            }
        }
    }

    fn update_low_qual_count(
        qual_map: &mut HashMap<usize, HashMap<DNA, usize>>,
        read: &SequenceRead,
        qual_threshold: u8,
    ) {
        fn add_low_qual_count(
            pos: usize,
            qual_count_map: &mut HashMap<usize, HashMap<DNA, usize>>,
            char: &DNA,
        ) {
            let inner_map = qual_count_map
                .entry(pos)
                .or_insert_with(DNA::create_a_dna_count_map);
            *inner_map.entry(char.clone()).or_insert(0) += 1;
        }
        if qual_threshold != 0_u8 {
            for (i, qual) in read.get_q_score_vec().iter().enumerate() {
                if *qual < qual_threshold {
                    let dna = read.get_seq_char(i);
                    add_low_qual_count(0, qual_map, &dna);
                    add_low_qual_count(i + 1, qual_map, &dna);
                }
            }
        }
    }

    /// update the map with a SequenceRead
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
}
