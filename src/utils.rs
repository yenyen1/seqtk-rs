use bio::io::{fasta, fastq};

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

