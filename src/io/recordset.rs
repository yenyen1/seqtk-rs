#[derive(Debug, Clone)]
pub struct RecordSetConfig {
    pub id_byte_limit: usize,
    pub seq_byte_limit: usize,
    pub qual_byte_limit: usize,
    pub record_size_limit: usize,
}
impl RecordSetConfig {
    pub fn fastq_default() -> Self {
        Self {
            id_byte_limit: 1024 * 1024,
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
    pub ids: Vec<u8>,
    pub seqs: Vec<u8>,
    pub quals: Vec<u8>,
    pub meta: Vec<(usize, usize, usize)>, // (id_len, seq_len, qual_len)
    config: RecordSetConfig,
}

impl OwnedRecordSet {
    pub fn new(config: RecordSetConfig) -> Self {
        Self {
            ids: Vec::with_capacity(config.id_byte_limit),
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
    pub fn push(&mut self, id: &[u8], seq: &[u8], qual: &[u8]) {
        self.ids.extend_from_slice(id);
        self.seqs.extend_from_slice(seq);
        let qual_len = if self.config.qual_byte_limit > 0 {
            self.quals.extend_from_slice(qual);
            qual.len()
        } else {
            0
        };
        self.meta.push((id.len(), seq.len(), qual_len));
    }
    pub fn clear(&mut self) {
        self.ids.clear();
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
            id_offset: 0,
            seq_offset: 0,
            qual_offset: 0,
        }
    }
}

pub struct RecordSlice<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

pub struct RecordSetIter<'a> {
    data: &'a OwnedRecordSet,
    index: usize,
    id_offset: usize,
    seq_offset: usize,
    qual_offset: usize,
}
impl<'a> Iterator for RecordSetIter<'a> {
    type Item = RecordSlice<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.data.meta.len() {
            return None;
        }

        let (id_len, seq_len, qual_len) = self.data.meta[self.index];
        let recode_slice = RecordSlice {
            id: &self.data.ids[self.id_offset..self.id_offset + id_len],
            seq: &self.data.seqs[self.seq_offset..self.seq_offset + seq_len],
            qual: &self.data.quals[self.qual_offset..self.qual_offset + qual_len],
        };

        self.index += 1;
        self.id_offset += id_len;
        self.seq_offset += seq_len;
        self.qual_offset += qual_len;

        Some(recode_slice)
    }
}
