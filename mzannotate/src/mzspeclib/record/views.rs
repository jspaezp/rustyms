use super::metadata::Ref;
use super::*;

#[derive(Clone, Copy, Debug)]
pub(super) enum Owner<'a> {
    Record(&'a SpectrumRecord<'a>),
    Header(LibraryContext<'a>),
}
impl<'a> Owner<'a> {
    fn raw(self, header: bool) -> &'a Raw {
        match self {
            Self::Record(r) if !header => &r.raw,
            Self::Record(r) => &r.context.metadata.raw,
            Self::Header(c) => &c.metadata.raw,
        }
    }
    fn identity(self) -> OwnerKey {
        match self {
            Self::Record(r) => OwnerKey::Record(r.owner, r.generation),
            Self::Header(c) => OwnerKey::Header(c.metadata.source),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OwnerKey {
    Record(RecordOwner, u64),
    Header(SourceId),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EffectiveGroupId {
    owner: OwnerKey,
    scope: usize,
    number: u64,
}
impl EffectiveGroupId {
    pub fn number(self) -> u64 {
        self.number
    }
}
#[derive(Clone, Copy, Debug)]
pub(super) enum GroupFilter {
    All,
    Group(Option<u64>),
}

#[derive(Clone, Copy, Debug)]
pub struct AttributeView<'a> {
    pub(super) owner: Owner<'a>,
    pub(super) scope: usize,
    pub(super) entries: Option<&'a [Ref]>,
    pub(super) group: GroupFilter,
}
impl<'a> AttributeView<'a> {
    pub(super) fn local(record: &'a SpectrumRecord<'a>, scope: usize) -> Self {
        Self {
            owner: Owner::Record(record),
            scope,
            entries: None,
            group: GroupFilter::All,
        }
    }
    pub fn iter(&self) -> AttributeIter<'a> {
        AttributeIter {
            view: *self,
            offset: 0,
        }
    }
    pub fn by_accession(&self, accession: Curie) -> impl Iterator<Item = AttributeOccurrence<'a>> {
        self.iter().filter(move |v| v.accession() == Some(accession))
    }
    pub fn group(&self, id: EffectiveGroupId) -> Option<Self> {
        if id.owner != self.owner.identity() || id.scope != self.scope {
            return None;
        }
        let view = Self {
            group: GroupFilter::Group(Some(id.number)),
            ..*self
        };
        (!view.is_empty()).then_some(view)
    }
    pub fn len(&self) -> usize {
        self.iter().count()
    }
    pub fn is_empty(&self) -> bool {
        self.iter().next().is_none()
    }
    /// Custom name/value groups, preserving duplicate names and all paired values.
    pub fn custom_values(&self, name: &'a str) -> impl Iterator<Item = CustomPropertyView<'a>> {
        self.by_accession(mzcv::curie!(MS:1003275))
            .filter(move |v| v.raw_value().trim_matches('"') == name)
            .map(|name| CustomPropertyView { name })
    }
}
impl<'a> IntoIterator for AttributeView<'a> {
    type Item = AttributeOccurrence<'a>;
    type IntoIter = AttributeIter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
#[derive(Debug)]
pub struct AttributeIter<'a> {
    view: AttributeView<'a>,
    offset: usize,
}
impl<'a> Iterator for AttributeIter<'a> {
    type Item = AttributeOccurrence<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let reference = if let Some(entries) = self.view.entries {
                *entries.get(self.offset)?
            } else {
                let raw = self.view.owner.raw(false);
                let scope = raw.scopes.get(self.view.scope)?;
                let index = scope.attrs.start + self.offset;
                if index >= scope.attrs.end {
                    return None;
                }
                Ref {
                    header: matches!(self.view.owner, Owner::Header(_)),
                    index,
                    group: raw.attrs[index].group.map(u64::from),
                    set: matches!(scope.id, ScopeId::AttributeSet(_)).then_some(self.view.scope),
                    claim: None,
                }
            };
            self.offset += 1;
            if matches!(self.view.group, GroupFilter::Group(g) if g != reference.group) {
                continue;
            }
            return Some(AttributeOccurrence {
                view: self.view,
                reference,
            });
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct AttributeOccurrence<'a> {
    view: AttributeView<'a>,
    reference: Ref,
}
impl<'a> AttributeOccurrence<'a> {
    fn raw(self) -> &'a Raw {
        self.view.owner.raw(self.reference.header)
    }
    fn span(self) -> &'a AttributeSpan {
        &self.raw().attrs[self.reference.index]
    }
    pub fn accession(&self) -> Option<Curie> {
        self.span().accession
    }
    pub fn name(&self) -> &'a str {
        &self.raw().text[self.span().name.clone()]
    }
    pub fn raw_value(&self) -> &'a str {
        &self.raw().text[self.span().value.clone()]
    }
    pub fn value(&self) -> Result<ValueView<'a>, &'a RecordError> {
        let raw = self.raw();
        let values = raw
            .values
            .get(|out| {
                for attr in &raw.attrs {
                    out.values.push(
                        StoredValue::parse(&raw.text, attr.value.clone(), &mut out.list_values)
                            .map_err(|message| {
                                RecordError::new(RecordErrorKind::Malformed, message, attr.position)
                            }),
                    );
                }
                Ok(())
            })
            .expect("value batch stores errors per occurrence");
        values.values[self.reference.index]
            .as_ref()
            .map(|v| v.view(&raw.text, &values.list_values))
    }
    pub fn to_f64(&self) -> Result<f64, RecordError> {
        match self.value() {
            Ok(ValueView::Integer(v)) => Ok(v as f64),
            Ok(ValueView::Float(v)) => Ok(v),
            Err(error) => Err(error.clone()),
            _ => Err(RecordError::new(
                RecordErrorKind::Malformed,
                "expected finite numeric attribute",
                self.span().position,
            )),
        }
    }
    pub fn group_id(&self) -> Option<EffectiveGroupId> {
        self.reference.group.map(|number| EffectiveGroupId {
            owner: self.view.owner.identity(),
            scope: self.view.scope,
            number,
        })
    }
    pub fn group(&self) -> Option<AttributeView<'a>> {
        self.group_id().and_then(|id| self.view.group(id))
    }
    pub fn origin(&self) -> AttributeOrigin<'a> {
        let header = self.view.owner.raw(true);
        AttributeOrigin {
            position: self.span().position,
            source_group: self.span().group,
            attribute_set: self
                .reference
                .set
                .and_then(|i| header.scopes[i].name.as_ref())
                .map(|r| &header.text[r.clone()]),
            claim_position: self
                .reference
                .claim
                .map(|i| self.view.owner.raw(false).attrs[i].position),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct AttributeOrigin<'a> {
    pub position: SourcePosition,
    pub source_group: Option<u32>,
    pub attribute_set: Option<&'a str>,
    pub claim_position: Option<SourcePosition>,
}

#[derive(Clone, Copy, Debug)]
pub enum ValueView<'a> {
    Empty,
    Integer(i64),
    Float(f64),
    Boolean(bool),
    Text(&'a str),
    Term(TermView<'a>),
    List(ListView<'a>),
}

#[derive(Debug, Default)]
pub(super) struct DecodedValues {
    values: Vec<Result<StoredValue, RecordError>>,
    list_values: Vec<StoredValue>,
}
impl Clear for DecodedValues {
    fn clear_reuse(&mut self) {
        self.values.clear();
        self.list_values.clear();
    }
}
#[derive(Debug)]
enum StoredValue {
    Empty,
    Integer(i64),
    Float(f64),
    Boolean(bool),
    Text(Range<usize>),
    Term(Option<Curie>, Range<usize>, Range<usize>),
    List(Range<usize>),
}
impl StoredValue {
    fn parse(raw: &str, range: Range<usize>, lists: &mut Vec<Self>) -> Result<Self, &'static str> {
        let original = &raw[range.clone()];
        let text = original.trim();
        let start = range.start + original.len() - original.trim_start().len();
        let range = start..start + text.len();
        if text.starts_with('"') {
            return if text.ends_with('"') && text.len() >= 2 {
                Ok(Self::Text(start + 1..range.end - 1))
            } else {
                Err("unclosed quoted value")
            };
        }
        if let Some((accession, name)) = text.split_once('|') {
            if accession.is_empty() || name.is_empty() {
                return Err("malformed term value");
            }
            return Ok(Self::Term(
                accession.parse().ok(),
                start..start + accession.len(),
                start + accession.len() + 1..range.end,
            ));
        }
        if text.contains(',')
            && text
                .split(',')
                .all(|s| s.trim().parse::<f64>().is_ok() || s.trim().parse::<bool>().is_ok())
        {
            let first = lists.len();
            let mut offset = start;
            for part in text.split(',') {
                let value = Self::parse(raw, offset..offset + part.len(), lists)?;
                lists.push(value);
                offset += part.len() + 1;
            }
            return Ok(Self::List(first..lists.len()));
        }
        Ok(if text.is_empty() {
            Self::Empty
        } else if let Ok(v) = text.parse::<i64>() {
            Self::Integer(v)
        } else if let Ok(v) = text.parse::<f64>() {
            if !v.is_finite() {
                return Err("nonfinite value");
            }
            Self::Float(v)
        } else if let Ok(v) = text.parse::<bool>() {
            Self::Boolean(v)
        } else {
            Self::Text(range)
        })
    }
    fn view<'a>(&self, raw: &'a str, lists: &'a [Self]) -> ValueView<'a> {
        match self {
            Self::Empty => ValueView::Empty,
            Self::Integer(v) => ValueView::Integer(*v),
            Self::Float(v) => ValueView::Float(*v),
            Self::Boolean(v) => ValueView::Boolean(*v),
            Self::Text(r) => ValueView::Text(&raw[r.clone()]),
            Self::Term(accession, id, name) => ValueView::Term(TermView {
                accession: *accession,
                raw_accession: &raw[id.clone()],
                name: &raw[name.clone()],
            }),
            Self::List(r) => ValueView::List(ListView {
                raw,
                values: &lists[r.clone()],
            }),
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct TermView<'a> {
    pub accession: Option<Curie>,
    pub raw_accession: &'a str,
    pub name: &'a str,
}
#[derive(Clone, Copy, Debug)]
pub struct ListView<'a> {
    raw: &'a str,
    values: &'a [StoredValue],
}
impl<'a> ListView<'a> {
    pub fn iter(&self) -> impl ExactSizeIterator<Item = ValueView<'a>> + '_ {
        self.values.iter().map(|v| v.view(self.raw, self.values))
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CustomPropertyView<'a> {
    name: AttributeOccurrence<'a>,
}
impl<'a> CustomPropertyView<'a> {
    pub fn name(&self) -> AttributeOccurrence<'a> {
        self.name
    }
    pub fn values(&self) -> impl Iterator<Item = AttributeOccurrence<'a>> {
        self.name
            .group()
            .into_iter()
            .flat_map(|g| g.into_iter())
            .filter(|v| v.accession() == Some(mzcv::curie!(MS:1003276)))
    }
    pub fn is_ambiguous(&self) -> bool {
        self.name.group().is_none_or(|g| {
            g.by_accession(mzcv::curie!(MS:1003275)).count() != 1
                || g.by_accession(mzcv::curie!(MS:1003276)).count() != 1
        })
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ScopeView<'a> {
    pub(super) record: &'a SpectrumRecord<'a>,
    pub(super) index: usize,
}
impl<'a> ScopeView<'a> {
    pub fn id(&self) -> ScopeId {
        self.record.raw.scopes[self.index].id
    }
    pub fn local_attributes(&self) -> AttributeView<'a> {
        AttributeView::local(self.record, self.index)
    }
    pub fn attributes(&self) -> Result<AttributeView<'a>, RecordError> {
        self.record.effective(self.index)
    }
    pub fn declaration(&self) -> &'a str {
        &self.record.raw.text[self.record.raw.scopes[self.index].declaration.clone()]
    }
    pub fn diagnostics(&self) -> Result<impl Iterator<Item = &'a RecordError>, RecordError> {
        Ok(self.record.resolve(self.index)?.diagnostics.iter())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct HeaderView<'a> {
    pub(super) context: LibraryContext<'a>,
}
impl<'a> HeaderView<'a> {
    pub fn attributes(&self) -> AttributeView<'a> {
        AttributeView {
            owner: Owner::Header(self.context),
            scope: 0,
            entries: None,
            group: GroupFilter::All,
        }
    }
    pub fn local_attributes(&self) -> AttributeView<'a> {
        self.attributes()
    }
    pub fn attribute_sets(&self) -> impl Iterator<Item = AttributeSetView<'a>> {
        let context = self.context;
        context
            .metadata
            .raw
            .scopes
            .iter()
            .enumerate()
            .filter(|(_, s)| s.name.is_some())
            .map(move |(index, _)| AttributeSetView { context, index })
    }
    pub fn source_path(&self) -> Option<&'a Path> {
        self.context.metadata.path.as_deref()
    }
}
#[derive(Clone, Copy, Debug)]
pub struct AttributeSetView<'a> {
    context: LibraryContext<'a>,
    index: usize,
}
impl<'a> AttributeSetView<'a> {
    pub fn name(&self) -> &'a str {
        &self.context.metadata.raw.text
            [self.context.metadata.raw.scopes[self.index].name.clone().unwrap()]
    }
    pub fn entry_type(&self) -> EntryType {
        self.context.metadata.raw.scopes[self.index].kind.unwrap()
    }
    pub fn attributes(&self) -> AttributeView<'a> {
        AttributeView {
            owner: Owner::Header(self.context),
            scope: self.index,
            entries: None,
            group: GroupFilter::All,
        }
    }
}
