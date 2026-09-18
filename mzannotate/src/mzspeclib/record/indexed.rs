use super::*;
use std::io::{Seek, SeekFrom};

#[derive(Debug)]
pub(super) struct IndexEntry {
    position: SourcePosition,
    key: Id,
    name: Option<Box<str>>,
    scan: Option<usize>,
}

impl<'a, R: BufRead + Seek> MzSpecLibRecordReader<'a, R> {
    /// Index literal spectrum names/scan numbers without decoding chemistry or peaks.
    /// Input must be a library beginning at byte zero. Reader position is preserved.
    pub fn build_index(&mut self) -> Result<usize, RecordError> {
        if let Some(index) = &self.input.index {
            return Ok(index.len());
        }
        let source = self.context.metadata.source;
        let at = SourcePosition {
            source,
            line: self.input.line,
            byte_offset: self.input.offset,
        };
        let saved = self
            .input
            .reader
            .stream_position()
            .map_err(|e| RecordError::new(RecordErrorKind::Io, e.to_string(), at))?;
        let result = (|| {
            self.input.reader.rewind()?;
            let mut entries: Vec<IndexEntry> = Vec::new();
            let mut text = String::new();
            let mut line = 0;
            let mut offset = 0;
            let mut spectrum_scope = false;
            loop {
                text.clear();
                let bytes = self.input.reader.read_line(&mut text)?;
                if bytes == 0 {
                    break;
                }
                let value = text.trim_end_matches(['\r', '\n']);
                if value.starts_with("<Spectrum=") {
                    spectrum_scope = true;
                    let position = SourcePosition {
                        source,
                        line,
                        byte_offset: offset,
                    };
                    if let Ok(key) = declaration(value, "<Spectrum=", position) {
                        entries.push(IndexEntry {
                            position,
                            key,
                            name: None,
                            scan: None,
                        });
                    } else {
                        spectrum_scope = false;
                    }
                } else if value.starts_with('<') {
                    spectrum_scope = false;
                } else if spectrum_scope {
                    if let Some(entry) = entries.last_mut() {
                        if let Some(name) = value.strip_prefix("MS:1003061|library spectrum name=")
                        {
                            entry.name = Some(name.into());
                        }
                        if let Some(scan) = value.strip_prefix("MS:1003057|scan number=") {
                            entry.scan = scan.trim().parse().ok();
                        }
                    }
                }
                line += 1;
                offset += bytes as u64;
            }
            Ok::<_, std::io::Error>(entries)
        })();
        self.input
            .reader
            .seek(SeekFrom::Start(saved))
            .map_err(|e| RecordError::new(RecordErrorKind::Io, e.to_string(), at))?;
        let index = result.map_err(|e| RecordError::new(RecordErrorKind::Io, e.to_string(), at))?;
        let len = index.len();
        self.input.index = Some(index);
        Ok(len)
    }
    fn read_position(
        &mut self,
        position: Option<SourcePosition>,
        record: &mut SpectrumRecord<'a>,
    ) -> Result<bool, RecordError> {
        record.reset(self.context);
        let Some(pos) = position else {
            return Ok(false);
        };
        self.input
            .reader
            .seek(SeekFrom::Start(pos.byte_offset))
            .map_err(|e| RecordError::new(RecordErrorKind::Io, e.to_string(), pos))?;
        self.input.pending.clear();
        self.input.line = pos.line;
        self.input.offset = pos.byte_offset;
        self.input.eof = false;
        self.read_into(record)
    }
    pub fn read_by_index_into(
        &mut self,
        index: usize,
        record: &mut SpectrumRecord<'a>,
    ) -> Result<bool, RecordError> {
        self.build_index()?;
        self.read_position(
            self.input.index.as_ref().unwrap().get(index).map(|v| v.position),
            record,
        )
    }
    pub fn read_by_key_into(
        &mut self,
        key: Id,
        record: &mut SpectrumRecord<'a>,
    ) -> Result<bool, RecordError> {
        self.build_index()?;
        self.read_position(
            self.input
                .index
                .as_ref()
                .unwrap()
                .iter()
                .rev()
                .find(|v| v.key == key)
                .map(|v| v.position),
            record,
        )
    }
    pub fn read_by_name_into(
        &mut self,
        name: &str,
        record: &mut SpectrumRecord<'a>,
    ) -> Result<bool, RecordError> {
        self.build_index()?;
        self.read_position(
            self.input
                .index
                .as_ref()
                .unwrap()
                .iter()
                .rev()
                .find(|v| v.name.as_deref() == Some(name))
                .map(|v| v.position),
            record,
        )
    }
    pub fn read_by_scan_number_into(
        &mut self,
        scan: usize,
        record: &mut SpectrumRecord<'a>,
    ) -> Result<bool, RecordError> {
        self.build_index()?;
        self.read_position(
            self.input
                .index
                .as_ref()
                .unwrap()
                .iter()
                .rev()
                .find(|v| v.scan == Some(scan))
                .map(|v| v.position),
            record,
        )
    }
}
