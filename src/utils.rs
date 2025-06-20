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
// pub enum FxWriter {
//     Fasta(fasta::Writer<Stdout>),
//     Fastq(fastq::Writer<Stdout>),
// }
// impl FxWriter {
//     pub fn new(is_fasta: bool) -> FxWriter {
//         let std_out = io::stdout();
//         if is_fasta {
//             FxWriter::Fasta(fasta::Writer::new(std_out))
//         } else {
//             FxWriter::Fastq(fastq::Writer::new(std_out))
//         }
//     }
//     pub fn write(
//         &mut self,
//         id: &str,
//         seq: &[u8],
//         desc: Option<&str>,
//         qual: &[u8],
//     ) -> io::Result<()> {
//         match self {
//             FxWriter::Fasta(ref mut out) => {
//                 out.write(id, desc, seq)?;
//             }
//             FxWriter::Fastq(ref mut out) => {
//                 out.write(id, desc, seq, qual)?;
//             }
//         }
//         Ok(())
//     }
// }
