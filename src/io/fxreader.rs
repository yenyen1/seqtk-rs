use crate::io::fxerr::FxError;
use crate::io::recordset::OwnedRecordSet;

use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind, Read};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;
use seq_io::fasta::{self, Record as _};
use seq_io::fastq::{self, Record as _};

type BoxedRead = Box<dyn Read + Send>;

/// FASTA/FASTQ Reader
pub enum FxReader {
    Fasta(fasta::Reader<BoxedRead>),
    Fastq(fastq::Reader<BoxedRead>),
}
impl FxReader {
    /// Constructs a new `FxReader` from the specified path
    ///
    /// The input file is automatically inspected to determine:
    /// - Whether the file is compressed (gzip supported)
    /// - Whether the format is FASTA or FASTQ
    ///
    /// This reader is built on top of the `seq_io` crate.
    /// It implements compression detection and uses an enum to wrapper
    /// to unify different reader types.
    ///
    /// It returns an appropriate reader accordingly
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The file cannot be opened (`NotFound`, `PermissionDenied`, etc.)
    /// - The file is empty (`UnexpectedEof`)
    /// - The compression format is not supported: bzip2, xz, zstd, etc. (`InvalidData`)
    /// - First non-empty line does not start with `b'>'` or `b'@'` (`InvalidData`)
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, FxError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        let buf = reader.fill_buf()?;
        if buf.starts_with(b"BZh") {
            return Err(FxError::Io(Error::new(
                ErrorKind::InvalidData,
                "Unsupported format: bzip2 not supported",
            )));
        } else if buf.starts_with(&[0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00]) {
            return Err(FxError::Io(Error::new(
                ErrorKind::InvalidData,
                "Unsupported format: xz not supported",
            )));
        } else if buf.starts_with(&[0x28, 0xB5, 0x2F, 0xFD]) {
            return Err(FxError::Io(Error::new(
                ErrorKind::InvalidData,
                "Unsupported format: zstd compression not supported",
            )));
        }
        let is_gzip = { buf.starts_with(&[0x1f, 0x8b]) };

        let decoder: BoxedRead = if is_gzip {
            Box::new(MultiGzDecoder::new(reader))
        } else {
            Box::new(reader)
        };

        let mut format_reader = BufReader::new(decoder);
        let first_byte = {
            let buf = format_reader.fill_buf()?;
            if buf.is_empty() {
                return Err(FxError::Io(Error::new(
                    ErrorKind::UnexpectedEof,
                    "File is empty",
                )));
            }
            buf[0]
        };

        let boxed_reader: BoxedRead = Box::new(format_reader);
        match first_byte {
            b'>' => Ok(FxReader::Fasta(fasta::Reader::new(boxed_reader))),
            b'@' => Ok(FxReader::Fastq(fastq::Reader::new(boxed_reader))),
            _ => Err(FxError::Io(Error::new(
                ErrorKind::InvalidData,
                "Unknown format: Unknown compressed format or File is neither FASTA nor FASTQ.",
            ))),
        }
    }
}

/// ## BatchReader
/// It is typically used in parallel processing to read data in batches,
/// which are then distributed to worker threads.
pub struct BatchReader {
    inner: FxReader,
    leftover_fa: Option<fasta::OwnedRecord>,
    leftover_fq: Option<fastq::OwnedRecord>,
}
impl BatchReader {
    /// Constructs a new `BatchReader` from an `FxReader`.
    ///
    /// ### Example
    /// ```
    /// let reader = FxReader::new(path)?;
    /// let batch_reader = BatchReader::new(reader);
    /// ```
    pub fn new(inner: FxReader) -> Self {
        Self {
            inner,
            leftover_fa: None,
            leftover_fq: None,
        }
    }

    /// Fill out the buffer (`OwnedRecordSet`) with records
    ///
    /// ### Returns
    /// - Ok(true) if buffer is not empty
    /// - Ok(false) if buffer is empty (EOF)
    /// - Err(e) if record is failed to parse by seq-io
    ///
    /// ### Example
    /// ```
    /// use crate::io::recordset::{OwnedRecordSet, RecordSetConfig};
    ///
    /// let reader = FxReader::new(path)?;
    /// let batch_reader = BatchReader::new(reader);
    /// let mut record_set = OwnedRecordSet::new(RecordSetConfig::fastq_default());
    /// match reader.fill_batch(&mut batch) {
    ///         Ok(true) => {
    ///                 ...
    ///         },
    ///         Ok(false) => break, // EOF
    ///         Err(e) => return Err(e), // Return Error
    ///         
    ///     }
    /// ```
    pub fn fill_batch(&mut self, batch: &mut OwnedRecordSet) -> Result<bool, FxError> {
        batch.clear();

        match &mut self.inner {
            FxReader::Fasta(reader) => {
                if let Some(record) = self.leftover_fa.take() {
                    batch.push(record.id_bytes(), record.seq(), &[]);
                }

                while let Some(record) = reader.next() {
                    match record {
                        Ok(rec) => {
                            let seq = rec.full_seq();
                            if batch.is_overload(seq.len()) {
                                self.leftover_fa = Some(rec.to_owned_record());
                                break;
                            }
                            batch.push(rec.id_bytes(), &seq, &[]);
                        }
                        Err(e) => {
                            // Error return from `seq-io`
                            return Err(FxError::from(e));
                        }
                    };
                }
            }

            FxReader::Fastq(reader) => {
                if let Some(record) = self.leftover_fq.take() {
                    batch.push(record.id_bytes(), record.seq(), &[]);
                }

                while let Some(record) = reader.next() {
                    match record {
                        Ok(rec) => {
                            if batch.is_overload(rec.seq().len()) {
                                self.leftover_fq = Some(rec.to_owned_record());
                            }
                            batch.push(rec.id_bytes(), rec.seq(), rec.qual());
                        }
                        Err(e) => {
                            // Error return from `seq-io`
                            return Err(FxError::from(e));
                        }
                    }
                }
            }
        }
        Ok(!batch.is_empty()) // return True if batch is not empty
    }
}
