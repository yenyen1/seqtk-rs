use std::io::BufRead;

pub struct FxIterator<R: BufRead> {
    reader: R,
    lines_per_chunk: usize,
}
impl<R: BufRead> FxIterator<R> {
    pub fn new(reader: R, lines_per_chunk: usize) -> Self {
        FxIterator {
            reader,
            lines_per_chunk,
        }
    }
}
impl<R: BufRead> Iterator for FxIterator<R> {
    type Item = Vec<String>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut lines = Vec::with_capacity(self.lines_per_chunk);

        for _ in 0..self.lines_per_chunk {
            let mut line = String::new();
            if self.reader.read_line(&mut line).unwrap_or(0) == 0 {
                if lines.is_empty() {
                    return None;
                } else {
                    return Some(lines);
                }
            }
            lines.push(line);
        }
        Some(lines)
    }
}
