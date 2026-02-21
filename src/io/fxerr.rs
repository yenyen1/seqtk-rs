use seq_io::fasta::Error as SeqioFaErr;
use seq_io::fastq::Error as SeqioFqErr;
use std::{fmt, io};

#[derive(Debug)]
pub enum FxError {
    Io(io::Error),
    SeqioFaErr(SeqioFaErr),
    SeqioFqErr(SeqioFqErr),
}

impl From<io::Error> for FxError {
    fn from(e: io::Error) -> Self {
        FxError::Io(e)
    }
}
impl From<SeqioFaErr> for FxError {
    fn from(e: SeqioFaErr) -> Self {
        match e {
            SeqioFaErr::Io(io_err) => FxError::Io(io_err),
            // ErrorTypes:
            // - InvalidStart (check b'>' at the beginning) should never occur (handled by `FxReader::new`)
            // - BufferLimit occurs when seqio buffer limitation reached
            other_err => FxError::SeqioFaErr(other_err),
        }
    }
}
impl From<SeqioFqErr> for FxError {
    fn from(e: SeqioFqErr) -> Self {
        match e {
            SeqioFqErr::Io(io_err) => FxError::Io(io_err),
            // ErrorType:
            // - BufferLimit occurs when seqio buffer limitation reached
            // - InvalidStart (not b'@'), UnequalLengths (seq.len != qual.len), InvalidSep (not '+') happend when calling `validate()` but not check every records
            // - UnexpectedEnd
            other_err => FxError::SeqioFqErr(other_err),
        }
    }
}
impl fmt::Display for FxError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FxError::Io(e) => e.fmt(f),
            FxError::SeqioFaErr(e) => e.fmt(f),
            FxError::SeqioFqErr(e) => e.fmt(f),
        }
    }
}
