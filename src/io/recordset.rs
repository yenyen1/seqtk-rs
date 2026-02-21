#[derive(Debug, Clone)]
pub struct RecordSetConfig {
    pub seq_byte_limit: usize,
    pub qual_byte_limit: usize,
    pub record_size_limit: usize,
}
impl RecordSetConfig {
    pub fn fastq_default() -> Self {
        Self {
            // id_byte_limit: 1024 * 1024,
            seq_byte_limit: 4 * 1024 * 1024,
            qual_byte_limit: 4 * 1024 * 1024,
            record_size_limit: 16 * 1024,
        }
    }
    pub fn fasta_default() -> Self {
        Self {
            qual_byte_limit: 0,
            ..Self::fastq_default()
        }
    }
}

/// ## OwnedRecordSet
/// It serves as a movable set of records that is typically used in parallel processing.
///
pub struct OwnedRecordSet {
    // pub ids: Vec<u8>,
    pub seqs: Vec<u8>,
    pub quals: Vec<u8>,
    pub meta: Vec<usize>, // (seq_len or qual_len)
    config: RecordSetConfig,
}

impl OwnedRecordSet {
    pub fn new(config: RecordSetConfig) -> Self {
        Self {
            seqs: Vec::with_capacity(config.seq_byte_limit),
            quals: Vec::with_capacity(config.qual_byte_limit),
            meta: Vec::with_capacity(config.record_size_limit),
            config,
        }
    }
    pub fn is_overload(&self, next_seq_len: usize) -> bool {
        if self.is_empty() {
            return false;
        }
        self.meta.len() >= self.config.record_size_limit
            || self.seqs.len() + next_seq_len > self.config.seq_byte_limit
    }
    pub fn push(&mut self, seq: &[u8], qual: &[u8]) {
        self.seqs.extend_from_slice(seq);
        if self.config.qual_byte_limit > 0 {
            self.quals.extend_from_slice(qual);
        }
        self.meta.push(seq.len());
    }
    pub fn clear(&mut self) {
        self.seqs.clear();
        self.quals.clear();
        self.meta.clear();
    }
    pub fn is_empty(&self) -> bool {
        self.meta.is_empty()
    }
    pub fn iter<'a>(&'a self) -> RecordSetIter<'a> {
        RecordSetIter {
            data: self,
            index: 0,
            seq_offset: 0,
            // qual_offset: 0,
        }
    }
    pub fn iter_seq<'a>(&'a self) -> SeqIter<'a> {
        SeqIter {
            data: self,
            index: 0,
            seq_offset: 0,
        }
    }
}

pub struct RecordSlice<'a> {
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

pub struct RecordSetIter<'a> {
    data: &'a OwnedRecordSet,
    index: usize,
    seq_offset: usize,
}
impl<'a> Iterator for RecordSetIter<'a> {
    type Item = RecordSlice<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.data.meta.len() {
            return None;
        }

        let seq_len = self.data.meta[self.index];
        let recode_slice = RecordSlice {
            seq: &self.data.seqs[self.seq_offset..self.seq_offset + seq_len],
            qual: &self.data.quals[self.seq_offset..self.seq_offset + seq_len],
        };

        self.index += 1;
        self.seq_offset += seq_len;

        Some(recode_slice)
    }
}

pub struct SeqSlice<'a> {
    pub seq: &'a [u8],
}

pub struct SeqIter<'a> {
    data: &'a OwnedRecordSet,
    index: usize,
    seq_offset: usize,
}
impl<'a> Iterator for SeqIter<'a> {
    type Item = SeqSlice<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.data.meta.len() {
            return None;
        }

        let seq_len = self.data.meta[self.index];
        let seq_slice = SeqSlice {
            seq: &self.data.seqs[self.seq_offset..self.seq_offset + seq_len],
        };

        self.index += 1;
        self.seq_offset += seq_len;

        Some(seq_slice)
    }
}
