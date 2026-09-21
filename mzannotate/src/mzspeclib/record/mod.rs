//! Owned, section-lazy mzSpecLib records with reusable range-backed storage.
//!
//! [`MzSpecLibLibrary`] owns the header and input. Records borrow its metadata,
//! own their raw spectrum text, and decode whole sections on first access.
//! Use [`MzSpecLibRecordReader::read_into`] to reuse allocations, or `records()`
//! for independent records. For worker-side structural parsing, use
//! [`MzSpecLibRecordReader::read_frame_into`] and [`SpectrumFrame::record`]. See `docs/mzspeclib-records.md` for usage and policies.
#![allow(missing_docs)]

mod decode;
mod indexed;
mod metadata;
#[cfg(test)]
mod tests;
mod views;
pub use decode::*;
pub use views::*;

use super::{EntryType, Id};
use mzcore::ontology::Ontologies;
use mzcv::Curie;
use std::{
    cell::{OnceCell, RefCell},
    fmt,
    io::BufRead,
    ops::Range,
    path::{Path, PathBuf},
    sync::Arc,
};

// Private constructors prevent consumers forging ownership tokens. These are
// process-local identities, not serialized keys, ordinals or content hashes.
mod identity {
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(1);
    fn next() -> u64 {
        NEXT.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |n| n.checked_add(1))
            .expect("ownership token space exhausted")
    }
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
    pub struct SourceId(u64);
    impl SourceId {
        pub(super) fn new() -> Self {
            Self(next())
        }
    }
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub(super) struct RecordOwner(u64);
    impl RecordOwner {
        pub(super) fn new() -> Self {
            Self(next())
        }
    }
}
use identity::RecordOwner;
pub use identity::SourceId;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SourcePosition {
    pub source: SourceId,
    pub line: u64,
    pub byte_offset: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecordErrorKind {
    Io,
    Structure,
    Malformed,
    Unsupported,
    Conflict,
    Unresolved,
    NoRecord,
}

#[derive(Debug)]
struct ErrorDetail {
    kind: RecordErrorKind,
    message: String,
    position: SourcePosition,
    path: Option<PathBuf>,
    source_line: Option<String>,
}

/// An owned diagnostic; cloning copies a shared handle, never source buffers.
#[derive(Clone, Debug)]
pub struct RecordError(Arc<ErrorDetail>);
impl RecordError {
    fn new(kind: RecordErrorKind, message: impl Into<String>, position: SourcePosition) -> Self {
        Self(Arc::new(ErrorDetail {
            kind,
            message: message.into(),
            position,
            path: None,
            source_line: None,
        }))
    }
    pub fn kind(&self) -> RecordErrorKind {
        self.0.kind
    }
    pub fn position(&self) -> SourcePosition {
        self.0.position
    }
    pub fn source_path(&self) -> Option<&Path> {
        self.0.path.as_deref()
    }
    pub fn source_line(&self) -> Option<&str> {
        self.0.source_line.as_deref()
    }
    fn with_context(mut self, path: Option<&Path>, text: &str, origin: u64) -> Self {
        if let Some(detail) = Arc::get_mut(&mut self.0) {
            detail.path = path.map(Path::to_path_buf);
            if let Some(offset) = detail
                .position
                .byte_offset
                .checked_sub(origin)
                .and_then(|n| usize::try_from(n).ok())
            {
                detail.source_line =
                    text.get(offset..).and_then(|s| s.lines().next()).map(str::to_owned);
            }
        }
        self
    }
}
impl fmt::Display for RecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (line {})", self.0.message, self.0.position.line + 1)
    }
}
impl std::error::Error for RecordError {}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ScopeId {
    Library,
    Spectrum(Id),
    Analyte(Id),
    Interpretation(Id),
    InterpretationMember(Id),
    AttributeSet(usize),
}

#[derive(Debug)]
struct AttributeSpan {
    group: Option<u32>,
    accession: Option<Curie>,
    name: Range<usize>,
    value: Range<usize>,
    position: SourcePosition,
}
#[derive(Debug)]
struct Scope {
    id: ScopeId,
    kind: Option<EntryType>,
    name: Option<Range<usize>>,
    attrs: Range<usize>,
    declaration: Range<usize>,
}

#[derive(Debug, Default)]
struct Raw {
    text: String,
    attrs: Vec<AttributeSpan>,
    scopes: Vec<Scope>,
    peaks: Range<usize>,
}
impl Raw {
    fn clear(&mut self) {
        self.text.clear();
        self.attrs.clear();
        self.scopes.clear();
        self.peaks = 0..0;
    }
}

trait Clear {
    fn clear_reuse(&mut self);
}
impl<T> Clear for Vec<T> {
    fn clear_reuse(&mut self) {
        self.clear();
    }
}

#[derive(Debug)]
struct Completed<T> {
    storage: T,
    outcome: Result<(), RecordError>,
}
#[derive(Debug)]
struct Reusable<T> {
    completed: OnceCell<Completed<T>>,
    spare: RefCell<Option<T>>,
}
impl<T: Default + Clear> Default for Reusable<T> {
    fn default() -> Self {
        Self {
            completed: OnceCell::new(),
            spare: RefCell::new(Some(T::default())),
        }
    }
}
impl<T: Default + Clear> Reusable<T> {
    fn get(&self, fill: impl FnOnce(&mut T) -> Result<(), RecordError>) -> Result<&T, RecordError> {
        let complete = self.completed.get_or_init(|| {
            let mut storage = self.spare.borrow_mut().take().unwrap_or_default();
            let outcome = fill(&mut storage);
            if outcome.is_err() {
                storage.clear_reuse();
            }
            Completed { storage, outcome }
        });
        complete.outcome.as_ref().map(|()| &complete.storage).map_err(Clone::clone)
    }
    fn reset(&mut self) {
        if let Some(mut complete) = self.completed.take() {
            complete.storage.clear_reuse();
            *self.spare.get_mut() = Some(complete.storage);
        } else if let Some(storage) = self.spare.get_mut() {
            storage.clear_reuse();
        }
    }
}

#[derive(Debug)]
struct LibraryMetadata {
    raw: Raw,
    values: DecodedValues,
    path: Option<PathBuf>,
    source: SourceId,
}
#[derive(Clone, Copy, Debug)]
struct LibraryContext<'a> {
    metadata: &'a LibraryMetadata,
    ontologies: &'a Ontologies,
}

#[derive(Debug)]
struct Input<R> {
    reader: R,
    pending: String,
    line: u64,
    offset: u64,
    pending_position: SourcePosition,
    eof: bool,
    index: Option<Vec<indexed::IndexEntry>>,
}
impl<R: BufRead> Input<R> {
    fn line(&mut self, buffer: &mut String) -> Result<Option<SourcePosition>, RecordError> {
        buffer.clear();
        if !self.pending.is_empty() {
            std::mem::swap(buffer, &mut self.pending);
            return Ok(Some(self.pending_position));
        }
        if self.eof {
            return Ok(None);
        }
        let pos = SourcePosition {
            source: self.pending_position.source,
            line: self.line,
            byte_offset: self.offset,
        };
        let count = self.reader.read_line(buffer).map_err(|e| {
            self.eof = true;
            RecordError::new(RecordErrorKind::Io, e.to_string(), pos)
        })?;
        if count == 0 {
            self.eof = true;
            return Ok(None);
        }
        self.line += 1;
        self.offset += count as u64;
        Ok(Some(pos))
    }
    fn put_back(&mut self, buffer: &mut String, position: SourcePosition) {
        std::mem::swap(&mut self.pending, buffer);
        self.pending_position = position;
    }
}

/// Library/header owner. The input is streamed; callers control retained record count.
/// Header values are decoded eagerly into immutable shared metadata. Malformed
/// values remain per-occurrence diagnostics rather than failing library opening.
#[derive(Debug)]
pub struct MzSpecLibLibrary<'o, R: BufRead> {
    metadata: LibraryMetadata,
    ontologies: &'o Ontologies,
    input: Input<R>,
}

impl<'o, R: BufRead> MzSpecLibLibrary<'o, R> {
    pub fn open(
        reader: R,
        path: Option<PathBuf>,
        ontologies: &'o Ontologies,
    ) -> Result<Self, RecordError> {
        let source = SourceId::new();
        let origin = SourcePosition {
            source,
            line: 0,
            byte_offset: 0,
        };
        let mut input = Input {
            reader,
            pending: String::new(),
            line: 0,
            offset: 0,
            pending_position: origin,
            eof: false,
            index: None,
        };
        let mut raw = Raw::default();
        let mut line = String::new();
        while let Some(position) = input.line(&mut line)? {
            if line.starts_with("<Spectrum=") {
                input.put_back(&mut line, position);
                break;
            }
            raw.text.push_str(&line);
        }
        scan(&mut raw, origin, true)?;
        let mut values = DecodedValues::default();
        values.decode(&raw);
        Ok(Self {
            metadata: LibraryMetadata {
                raw,
                values,
                path,
                source,
            },
            ontologies,
            input,
        })
    }
    pub fn header(&self) -> HeaderView<'_> {
        HeaderView {
            context: LibraryContext {
                metadata: &self.metadata,
                ontologies: self.ontologies,
            },
        }
    }
    pub fn reader(&mut self) -> MzSpecLibRecordReader<'_, R> {
        MzSpecLibRecordReader {
            input: &mut self.input,
            context: LibraryContext {
                metadata: &self.metadata,
                ontologies: self.ontologies,
            },
            line: String::new(),
        }
    }
}

#[derive(Debug)]
pub struct MzSpecLibRecordReader<'a, R: BufRead> {
    input: &'a mut Input<R>,
    context: LibraryContext<'a>,
    line: String,
}
impl<'a, R: BufRead> MzSpecLibRecordReader<'a, R> {
    pub fn header(&self) -> HeaderView<'a> {
        HeaderView {
            context: self.context,
        }
    }
    pub fn empty_record(&self) -> SpectrumRecord<'a> {
        SpectrumRecord::empty(self.context)
    }
    /// Allocate an empty reusable frame without parsing a spectrum.
    pub fn empty_frame(&self) -> SpectrumFrame<'a> {
        SpectrumFrame {
            record: self.empty_record(),
            parsed: None,
        }
    }
    /// Buffer a complete spectrum without parsing IDs, attributes or peak fields.
    /// Returns false at EOF. Structural errors are deferred to `frame.record()`;
    /// IO/UTF-8 errors remain reader errors. Existing borrowed views prevent refill.
    pub fn read_frame_into(&mut self, frame: &mut SpectrumFrame<'a>) -> Result<bool, RecordError> {
        frame.parsed = None;
        frame.record.reset_frame(self.context);
        let outcome = self.load_frame(&mut frame.record).map_err(|e| {
            e.with_context(
                self.context.metadata.path.as_deref(),
                &frame.record.raw.text,
                frame.record.origin.byte_offset,
            )
        });
        if outcome.is_err() {
            frame.record.reset_frame(self.context);
        }
        outcome
    }
    pub fn records(
        &mut self,
    ) -> impl Iterator<Item = Result<SpectrumRecord<'a>, RecordError>> + '_ {
        std::iter::from_fn(move || {
            let mut record = self.empty_record();
            match self.read_into(&mut record) {
                Ok(true) => Some(Ok(record)),
                Ok(false) => None,
                Err(e) => Some(Err(e)),
            }
        })
    }
    pub fn read_into(&mut self, record: &mut SpectrumRecord<'a>) -> Result<bool, RecordError> {
        record.reset(self.context);
        let outcome = self
            .load_frame(record)
            .and_then(|loaded| {
                if loaded {
                    record.parse_structure()?;
                }
                Ok(loaded)
            })
            .map_err(|e| {
                e.with_context(
                    self.context.metadata.path.as_deref(),
                    &record.raw.text,
                    record.origin.byte_offset,
                )
            });
        if outcome.is_err() {
            record.reset(self.context);
        }
        outcome
    }
    fn load_frame(&mut self, record: &mut SpectrumRecord<'a>) -> Result<bool, RecordError> {
        // On failure the previous record was fully consumed. Synchronize to the next declaration.
        loop {
            let Some(position) = self.input.line(&mut self.line)? else {
                return Ok(false);
            };
            if self.line.starts_with("<Spectrum=") {
                record.origin = position;
                record.raw.text.push_str(&self.line);
                break;
            }
        }
        while let Some(position) = self.input.line(&mut self.line)? {
            if self.line.starts_with("<Spectrum=") {
                self.input.put_back(&mut self.line, position);
                break;
            }
            record.raw.text.push_str(&self.line);
        }
        record.loaded = true;
        Ok(true)
    }
}

/// Reusable unparsed spectrum text plus retained storage for consumer-side decoding.
///
/// The reader only frames text. `record()` indexes it once, on the calling thread,
/// without copying its raw buffer. Later calls reuse the index or cached error.
/// All decoded views borrow this frame and prevent it from being refilled.
///
/// ```compile_fail
/// use mzannotate::mzspeclib::record::{MzSpecLibRecordReader, SpectrumFrame};
/// use std::io::Cursor;
/// fn invalid<'a>(reader: &mut MzSpecLibRecordReader<'a, Cursor<Vec<u8>>>,
///                frame: &mut SpectrumFrame<'a>) {
///     let record = frame.record().unwrap();
///     reader.read_frame_into(frame).unwrap();
///     println!("{:?}", record.key());
/// }
/// ```
#[derive(Debug)]
pub struct SpectrumFrame<'a> {
    record: SpectrumRecord<'a>,
    parsed: Option<Result<(), RecordError>>,
}
impl<'a> SpectrumFrame<'a> {
    pub fn raw_text(&self) -> &str {
        self.record.raw_text()
    }
    pub fn source_position(&self) -> Option<SourcePosition> {
        self.record.source_position()
    }
    pub fn source_path(&self) -> Option<&Path> {
        self.record.source_path()
    }
    pub fn header(&self) -> HeaderView<'_> {
        self.record.header()
    }
    /// Index structural metadata on this thread. Chemistry and peaks stay lazy.
    /// Structural failures preserve the raw text and are cached until refill.
    pub fn record(&mut self) -> Result<&SpectrumRecord<'a>, RecordError> {
        self.record.require()?;
        let outcome = self.parsed.get_or_insert_with(|| {
            self.record.reset_decoded();
            self.record.parse_structure().map_err(|e| {
                e.with_context(
                    self.record.context.metadata.path.as_deref(),
                    &self.record.raw.text,
                    self.record.origin.byte_offset,
                )
            })
        });
        outcome.as_ref().map_err(Clone::clone)?;
        Ok(&self.record)
    }
}

/// Owns spectrum text and decoded caches; references only immutable library metadata.
///
/// Live views prevent buffer reuse:
/// ```compile_fail
/// use mzannotate::mzspeclib::record::{MzSpecLibRecordReader, SpectrumRecord};
/// use std::io::Cursor;
/// fn invalid<'a>(reader: &mut MzSpecLibRecordReader<'a, Cursor<Vec<u8>>>,
///                record: &mut SpectrumRecord<'a>) {
///     let peak = record.peaks().unwrap().get(0).unwrap();
///     reader.read_into(record).unwrap();
///     println!("{}", peak.mz());
/// }
/// ```
/// Records also cannot outlive their header owner:
/// ```compile_fail
/// use mzannotate::mzspeclib::record::MzSpecLibLibrary;
/// use mzcore::ontology::STATIC_ONTOLOGIES;
/// let record = {
///     let mut library = MzSpecLibLibrary::open(
///         &b"<mzSpecLib>\n<Spectrum=1>\n<Peaks>\n"[..], None,
///         &STATIC_ONTOLOGIES).unwrap();
///     let mut reader = library.reader();
///     reader.records().next().unwrap().unwrap()
/// };
/// println!("{:?}", record.key());
/// ```
#[derive(Debug)]
pub struct SpectrumRecord<'a> {
    context: LibraryContext<'a>,
    raw: Raw,
    values: Reusable<DecodedValues>,
    origin: SourcePosition,
    owner: RecordOwner,
    generation: u64,
    loaded: bool,
    metadata: Vec<Reusable<metadata::ResolvedScope>>,
    analytes: Reusable<DecodedAnalytes>,
    interpretations: Reusable<DecodedInterpretations>,
    peaks: Reusable<DecodedPeaks>,
    annotations: Reusable<AnnotationStorage>,
}
impl<'a> SpectrumRecord<'a> {
    fn parse_structure(&mut self) -> Result<(), RecordError> {
        scan(&mut self.raw, self.origin, false)?;
        while self.metadata.len() < self.raw.scopes.len() {
            self.metadata.push(Reusable::default());
        }
        Ok(())
    }
    fn error(
        &self,
        kind: RecordErrorKind,
        message: impl Into<String>,
        position: SourcePosition,
    ) -> RecordError {
        let (text, offset) = if position.byte_offset < self.origin.byte_offset {
            (&self.context.metadata.raw.text, 0)
        } else {
            (&self.raw.text, self.origin.byte_offset)
        };
        RecordError::new(kind, message, position).with_context(
            self.context.metadata.path.as_deref(),
            text,
            offset,
        )
    }
    fn empty(context: LibraryContext<'a>) -> Self {
        Self {
            context,
            raw: Raw::default(),
            values: Reusable::default(),
            origin: SourcePosition {
                source: context.metadata.source,
                line: 0,
                byte_offset: 0,
            },
            owner: RecordOwner::new(),
            generation: 0,
            loaded: false,
            metadata: Vec::new(),
            analytes: Reusable::default(),
            interpretations: Reusable::default(),
            peaks: Reusable::default(),
            annotations: Reusable::default(),
        }
    }
    fn reset(&mut self, context: LibraryContext<'a>) {
        self.reset_decoded();
        self.reset_frame(context);
    }
    fn reset_decoded(&mut self) {
        for slot in &mut self.metadata {
            slot.reset();
        }
        self.analytes.reset();
        self.interpretations.reset();
        self.peaks.reset();
        self.annotations.reset();
        self.values.reset();
    }
    fn reset_frame(&mut self, context: LibraryContext<'a>) {
        self.raw.clear();
        self.context = context;
        self.loaded = false;
        if let Some(next) = self.generation.checked_add(1) {
            self.generation = next;
        } else {
            self.owner = RecordOwner::new();
            self.generation = 0;
        }
    }
    fn require(&self) -> Result<(), RecordError> {
        if self.loaded {
            Ok(())
        } else {
            Err(RecordError::new(
                RecordErrorKind::NoRecord,
                "no loaded spectrum",
                self.origin,
            ))
        }
    }
    pub fn is_loaded(&self) -> bool {
        self.loaded
    }
    pub fn key(&self) -> Option<Id> {
        self.raw.scopes.first().and_then(|s| match s.id {
            ScopeId::Spectrum(id) => Some(id),
            _ => None,
        })
    }
    pub fn source_position(&self) -> Option<SourcePosition> {
        self.loaded.then_some(self.origin)
    }
    pub fn source_path(&self) -> Option<&Path> {
        self.context.metadata.path.as_deref()
    }
    pub fn raw_text(&self) -> &str {
        &self.raw.text
    }
    pub fn header(&self) -> HeaderView<'_> {
        HeaderView {
            context: self.context,
        }
    }
    pub fn local_attributes(&self) -> AttributeView<'_> {
        AttributeView::local(self, 0)
    }
    pub fn attributes(&self) -> Result<AttributeView<'_>, RecordError> {
        self.require()?;
        self.effective(0)
    }
    pub fn analyte_scopes(&self) -> impl Iterator<Item = ScopeView<'_>> {
        self.raw
            .scopes
            .iter()
            .enumerate()
            .filter(|(_, s)| matches!(s.id, ScopeId::Analyte(_)))
            .map(|(index, _)| ScopeView {
                record: self,
                index,
            })
    }
    pub fn interpretation_scopes(&self) -> impl Iterator<Item = ScopeView<'_>> {
        self.raw
            .scopes
            .iter()
            .enumerate()
            .filter(|(_, s)| matches!(s.id, ScopeId::Interpretation(_)))
            .map(|(index, _)| ScopeView {
                record: self,
                index,
            })
    }
    pub fn scopes(&self) -> impl Iterator<Item = ScopeView<'_>> {
        (0..self.raw.scopes.len()).map(|index| ScopeView {
            record: self,
            index,
        })
    }
    pub fn attribute_set_claims(&self) -> impl Iterator<Item = AttributeOccurrence<'_>> {
        self.local_attributes()
            .into_iter()
            .filter(|a| a.accession() == Some(mzcv::curie!(MS:1003212)))
    }
    pub fn resolution_diagnostics(
        &self,
    ) -> Result<impl Iterator<Item = &RecordError>, RecordError> {
        self.require()?;
        Ok(self.resolve(0)?.diagnostics.iter())
    }
}

fn declaration(line: &str, prefix: &str, pos: SourcePosition) -> Result<Id, RecordError> {
    line.strip_prefix(prefix)
        .and_then(|s| s.strip_suffix('>'))
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| RecordError::new(RecordErrorKind::Structure, "invalid declaration", pos))
}

fn scan(raw: &mut Raw, origin: SourcePosition, header: bool) -> Result<(), RecordError> {
    let mut offset = 0;
    let mut peaks = false;
    for (line_number, full) in raw.text.split_inclusive('\n').enumerate() {
        let line = full.trim_end_matches(['\r', '\n']);
        let position = SourcePosition {
            source: origin.source,
            line: origin.line + line_number as u64,
            byte_offset: origin.byte_offset + offset as u64,
        };
        let end = offset + full.len();
        if line.trim().is_empty() || line.starts_with('#') {
            offset = end;
            continue;
        }
        if line == "<Peaks>" && !header {
            peaks = true;
            raw.peaks = end..raw.text.len();
            offset = end;
            continue;
        }
        if peaks {
            offset = end;
            continue;
        }
        if line.starts_with('<') {
            let (id, kind, name) = if header && line == "<mzSpecLib>" {
                (ScopeId::Library, None, None)
            } else if header && line.starts_with("<AttributeSet ") {
                let body =
                    line.strip_prefix("<AttributeSet ").unwrap().strip_suffix('>').ok_or_else(
                        || {
                            RecordError::new(
                                RecordErrorKind::Structure,
                                "invalid attribute set",
                                position,
                            )
                        },
                    )?;
                let (kind, name) = body.split_once('=').ok_or_else(|| {
                    RecordError::new(RecordErrorKind::Structure, "missing set name", position)
                })?;
                let kind = match kind {
                    "Spectrum" => EntryType::Spectrum,
                    "Analyte" => EntryType::Analyte,
                    "Interpretation" => EntryType::Interpretation,
                    "Cluster" => EntryType::Cluster,
                    _ => {
                        return Err(RecordError::new(
                            RecordErrorKind::Unsupported,
                            "unknown attribute set scope",
                            position,
                        ));
                    }
                };
                (
                    ScopeId::AttributeSet(raw.scopes.len()),
                    Some(kind),
                    Some(offset + line.len() - 1 - name.len()..offset + line.len() - 1),
                )
            } else if !header && line.starts_with("<Spectrum=") {
                (
                    ScopeId::Spectrum(declaration(line, "<Spectrum=", position)?),
                    Some(EntryType::Spectrum),
                    None,
                )
            } else if !header && line.starts_with("<Analyte=") {
                let id = declaration(line, "<Analyte=", position)?;
                if id == 0 {
                    return Err(RecordError::new(
                        RecordErrorKind::Structure,
                        "analyte ID cannot be zero",
                        position,
                    ));
                }
                (ScopeId::Analyte(id), Some(EntryType::Analyte), None)
            } else if !header && line.starts_with("<Interpretation=") {
                (
                    ScopeId::Interpretation(declaration(line, "<Interpretation=", position)?),
                    Some(EntryType::Interpretation),
                    None,
                )
            } else if !header && line.starts_with("<InterpretationMember=") {
                (
                    ScopeId::InterpretationMember(declaration(
                        line,
                        "<InterpretationMember=",
                        position,
                    )?),
                    None,
                    None,
                )
            } else {
                return Err(RecordError::new(
                    RecordErrorKind::Structure,
                    "unexpected declaration",
                    position,
                ));
            };
            raw.scopes.push(Scope {
                id,
                kind,
                name,
                attrs: raw.attrs.len()..raw.attrs.len(),
                declaration: offset..end,
            });
        } else {
            let scope = raw.scopes.last_mut().ok_or_else(|| {
                RecordError::new(
                    RecordErrorKind::Structure,
                    "attribute before declaration",
                    position,
                )
            })?;
            let (group, start) = if let Some(rest) = line.strip_prefix('[') {
                let (id, _) = rest.split_once(']').ok_or_else(|| {
                    RecordError::new(
                        RecordErrorKind::Structure,
                        "unclosed attribute group",
                        position,
                    )
                })?;
                (
                    Some(id.parse().map_err(|_| {
                        RecordError::new(RecordErrorKind::Structure, "invalid group ID", position)
                    })?),
                    id.len() + 2,
                )
            } else {
                (None, 0)
            };
            let equals = line[start..].find('=').map(|n| start + n).ok_or_else(|| {
                RecordError::new(
                    RecordErrorKind::Structure,
                    "missing attribute equals",
                    position,
                )
            })?;
            let name = &line[start..equals];
            let accession = name.split('|').next().and_then(|s| s.parse().ok());
            raw.attrs.push(AttributeSpan {
                group,
                accession,
                name: offset + start..offset + equals,
                value: offset + equals + 1..offset + line.len(),
                position,
            });
            scope.attrs.end = raw.attrs.len();
        }
        offset = end;
    }
    if raw.scopes.is_empty() || (header && raw.scopes[0].id != ScopeId::Library) {
        return Err(RecordError::new(
            RecordErrorKind::Structure,
            "missing library/record declaration",
            origin,
        ));
    }
    if !header && !peaks {
        return Err(RecordError::new(
            RecordErrorKind::Structure,
            "missing Peaks section",
            origin,
        ));
    }
    Ok(())
}
