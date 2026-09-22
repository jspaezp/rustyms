use super::*;
use mzcv::curie;

#[derive(Clone, Copy, Debug)]
pub(super) struct Ref {
    pub header: bool,
    pub index: usize,
    pub group: Option<u64>,
    pub set: Option<usize>,
    pub claim: Option<usize>,
}

#[derive(Debug, Default)]
pub(super) struct ResolvedScope {
    pub entries: Vec<Ref>,
    pub diagnostics: Vec<RecordError>,
}
impl Clear for ResolvedScope {
    fn clear_reuse(&mut self) {
        self.entries.clear();
        self.diagnostics.clear();
    }
}

impl SpectrumRecord<'_> {
    pub(super) fn resolve(&self, index: usize) -> Result<&ResolvedScope, RecordError> {
        self.metadata[index].get(|out| {
            let scope = &self.raw.scopes[index];
            let Some(kind) = scope.kind else {
                if matches!(scope.id, ScopeId::InterpretationMember(_)) {
                    out.diagnostics.push(self.error(
                        RecordErrorKind::Unsupported,
                        "interpretation member retained as raw evidence",
                        self.origin,
                    ));
                }
                out.entries.extend(scope.attrs.clone().map(|index| Ref {
                    header: false,
                    index,
                    group: self.raw.attrs[index].group.map(u64::from),
                    set: None,
                    claim: None,
                }));
                return Ok(());
            };
            let mut inherited_group = scope
                .attrs
                .clone()
                .filter_map(|i| self.raw.attrs[i].group)
                .max()
                .map_or(0, |g| u64::from(g) + 1);
            // Ungrouped default, named claims in claim order, then local values.
            self.apply_set(out, kind, "all", None, None, &mut inherited_group, false);
            for i in scope.attrs.clone().filter(|i| self.raw.attrs[*i].group.is_none()) {
                if self.raw.attrs[i].accession == Some(curie!(MS:1003212)) {
                    self.apply_set(
                        out,
                        kind,
                        self.raw.text[self.raw.attrs[i].value.clone()].trim(),
                        Some(i),
                        None,
                        &mut inherited_group,
                        true,
                    );
                }
            }
            self.merge_source(out, false, index, None, None, None, None);
            for i in scope.attrs.clone() {
                let Some(group) = self.raw.attrs[i].group else {
                    continue;
                };
                if scope
                    .attrs
                    .clone()
                    .take_while(|j| *j < i)
                    .any(|j| self.raw.attrs[j].group == Some(group))
                {
                    continue;
                }
                for claim in scope.attrs.clone().filter(|j| {
                    self.raw.attrs[*j].group == Some(group)
                        && self.raw.attrs[*j].accession == Some(curie!(MS:1003212))
                }) {
                    self.apply_set(
                        out,
                        kind,
                        self.raw.text[self.raw.attrs[claim].value.clone()].trim(),
                        Some(claim),
                        Some(u64::from(group)),
                        &mut inherited_group,
                        true,
                    );
                }
                self.merge_source(
                    out,
                    false,
                    index,
                    Some(group),
                    Some(u64::from(group)),
                    None,
                    None,
                );
            }
            Ok(())
        })
    }
    pub(super) fn effective(&self, scope: usize) -> Result<AttributeView<'_>, RecordError> {
        let entries = &self.resolve(scope)?.entries;
        Ok(AttributeView {
            owner: Owner::Record(self),
            scope,
            entries: Some(entries),
            group: GroupFilter::All,
        })
    }
    fn apply_set(
        &self,
        out: &mut ResolvedScope,
        kind: EntryType,
        name: &str,
        claim: Option<usize>,
        group: Option<u64>,
        next_group: &mut u64,
        diagnose: bool,
    ) {
        let header = &self.context.metadata.raw;
        let mut found = false;
        for (set, scope) in header.scopes.iter().enumerate().filter(|(_, s)| {
            s.kind == Some(kind) && s.name.as_ref().is_some_and(|r| &header.text[r.clone()] == name)
        }) {
            found = true;
            self.merge_source(out, true, set, None, group, Some(set), claim);
            // A grouped claim imports only the set's ungrouped properties.
            if group.is_none() {
                for index in scope.attrs.clone() {
                    let Some(source_group) = header.attrs[index].group else {
                        continue;
                    };
                    if scope
                        .attrs
                        .clone()
                        .take_while(|j| *j < index)
                        .any(|j| header.attrs[j].group == Some(source_group))
                    {
                        continue;
                    }
                    self.merge_source(
                        out,
                        true,
                        set,
                        Some(source_group),
                        Some(*next_group),
                        Some(set),
                        claim,
                    );
                    *next_group += 1;
                }
            }
        }
        if !found && diagnose {
            out.diagnostics.push(self.error(
                RecordErrorKind::Unresolved,
                format!("undefined attribute set {name}"),
                claim.map_or(self.origin, |i| self.raw.attrs[i].position),
            ));
        }
    }
    fn merge_source(
        &self,
        out: &mut ResolvedScope,
        header: bool,
        scope: usize,
        source_group: Option<u32>,
        group: Option<u64>,
        set: Option<usize>,
        claim: Option<usize>,
    ) {
        let raw = if header {
            &self.context.metadata.raw
        } else {
            &self.raw
        };
        let attrs = raw.scopes[scope].attrs.clone();
        for index in attrs.clone().filter(|i| raw.attrs[*i].group == source_group) {
            // Remove lower-priority occurrences together, never same-source repeats.
            let first = !attrs
                .clone()
                .take_while(|i| *i < index)
                .any(|i| raw.attrs[i].group == source_group && same_key(raw, i, raw, index));
            let matches = |old: &Ref| old.group == group && same_key(
                if old.header { &self.context.metadata.raw } else { &self.raw }, old.index, raw, index);
            let insertion = if first {
                let insertion = out.entries.iter().position(&matches).unwrap_or(out.entries.len());
                out.entries.retain(|old| !matches(old));
                insertion
            } else {
                out.entries.iter().rposition(matches).map_or(out.entries.len(), |i| i + 1)
            };
            out.entries.insert(insertion, Ref { header, index, group, set, claim });
        }
    }
}

fn same_key(a: &Raw, i: usize, b: &Raw, j: usize) -> bool {
    match (a.attrs[i].accession, b.attrs[j].accession) {
        (Some(x), Some(y)) => x == y,
        _ => a.text[a.attrs[i].name.clone()] == b.text[b.attrs[j].name.clone()],
    }
}
