use super::*;
use crate::mzspeclib::{
    Analyte, AnalyteTarget, Attribute, Attributes, Interpretation, ProteinDescription,
    populate_spectrum_description_from_attributes,
};
use crate::{
    fragment::{Fragment, PeakAnnotation},
    spectrum::{AnnotatedPeak, AnnotatedSpectrum},
};
use context_error::Context;
use mzcore::{
    chemistry::OutputMolecularFormula,
    prelude::*,
    system::{MassOverCharge, isize::Charge},
};
use mzcv::curie;
use mzdata::{
    mzpeaks::prelude::PeakCollectionMut,
    spectrum::{Precursor, ScanEvent, ScanWindow, SelectedIon},
};

/// Decoded analyte chemistry. Supplied protein and custom properties remain
/// borrowed through the corresponding analyte scope; owned export is explicit.
pub struct DecodedAnalyte {
    pub id: std::num::NonZeroU32,
    pub target: AnalyteTarget,
    scratch: mzcore::sequence::ProFormaScratch,
    formulas: Reusable<mzcore::sequence::FormulaBuffer>,
    spare_ion: PeptidoformIon,
    spare_formula: MolecularFormula,
}
impl std::fmt::Debug for DecodedAnalyte {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DecodedAnalyte")
            .field("id", &self.id)
            .field("target", &self.target)
            .finish()
    }
}
impl Clear for mzcore::sequence::FormulaBuffer {
    fn clear_reuse(&mut self) {
        self.clear();
    }
}
#[derive(Clone, Copy)]
enum TargetKind {
    Unknown,
    Peptidoform,
    Formula,
}
impl DecodedAnalyte {
    fn select_target(&mut self, kind: TargetKind) {
        if matches!(
            (&self.target, kind),
            (AnalyteTarget::PeptidoformIon(_), TargetKind::Peptidoform)
                | (AnalyteTarget::MolecularFormula(_), TargetKind::Formula)
        ) {
            return;
        }
        match std::mem::take(&mut self.target) {
            AnalyteTarget::PeptidoformIon(ion) => self.spare_ion = ion,
            AnalyteTarget::MolecularFormula(formula) => self.spare_formula = formula,
            AnalyteTarget::Unknown(_) => (),
        }
        self.target = match kind {
            TargetKind::Peptidoform => {
                AnalyteTarget::PeptidoformIon(std::mem::take(&mut self.spare_ion))
            }
            TargetKind::Formula => {
                AnalyteTarget::MolecularFormula(std::mem::take(&mut self.spare_formula))
            }
            TargetKind::Unknown => AnalyteTarget::Unknown(None),
        };
    }

    /// Calculate and borrow all formula alternatives using retained output buffers.
    /// General/ambiguous chemistry may still allocate inside the resolver.
    pub fn formulas(&self) -> &[MolecularFormula] {
        self.formulas
            .get(|out| {
                match &self.target {
                    AnalyteTarget::PeptidoformIon(ion) => {
                        out.calculate(ion);
                    }
                    AnalyteTarget::MolecularFormula(formula) => out.set_formula(formula),
                    AnalyteTarget::Unknown(_) => out.clear(),
                }
                Ok(())
            })
            .expect("formula cache initialization is infallible")
            .as_slice()
    }
}

#[derive(Debug, Default)]
pub(super) struct DecodedAnalytes {
    slots: Vec<DecodedAnalyte>,
    active: usize,
}
impl Clear for DecodedAnalytes {
    fn clear_reuse(&mut self) {
        self.active = 0;
        for slot in &mut self.slots {
            slot.formulas.reset();
        }
    }
}
impl DecodedAnalytes {
    fn as_slice(&self) -> &[DecodedAnalyte] {
        &self.slots[..self.active]
    }
    fn into_active(mut self) -> impl Iterator<Item = DecodedAnalyte> {
        self.slots.truncate(self.active);
        self.slots.into_iter()
    }
}

/// Decoded interpretation links. Other properties remain in borrowed scope views.
#[derive(Debug, Default)]
pub struct DecodedInterpretation {
    pub id: u32,
    pub probability: Option<f64>,
    pub analyte_refs: Vec<u32>,
    scope_index: usize,
}
#[derive(Debug, Default)]
pub(super) struct DecodedInterpretations {
    slots: Vec<DecodedInterpretation>,
    active: usize,
}
impl Clear for DecodedInterpretations {
    fn clear_reuse(&mut self) {
        self.active = 0;
    }
}
impl DecodedInterpretations {
    fn as_slice(&self) -> &[DecodedInterpretation] {
        &self.slots[..self.active]
    }
    fn into_active(mut self) -> impl Iterator<Item = DecodedInterpretation> {
        self.slots.truncate(self.active);
        self.slots.into_iter()
    }
}

#[derive(Debug, Default)]
pub(super) struct DecodedPeaks {
    mz: Vec<f64>,
    intensity: Vec<f32>,
    rows: Vec<Row>,
}
#[derive(Debug)]
struct Row {
    annotation: Option<Range<usize>>,
    aggregation: Option<Range<usize>>,
    extras: Option<Range<usize>>,
    position: SourcePosition,
}
impl Clear for DecodedPeaks {
    fn clear_reuse(&mut self) {
        self.mz.clear();
        self.intensity.clear();
        self.rows.clear();
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SuppliedField<'a> {
    Absent,
    Present(&'a str),
}
#[derive(Clone, Copy, Debug)]
pub struct PeaksView<'a> {
    record: &'a SpectrumRecord<'a>,
    storage: &'a DecodedPeaks,
}
impl<'a> PeaksView<'a> {
    pub fn len(&self) -> usize {
        self.storage.rows.len()
    }
    pub fn is_empty(&self) -> bool {
        self.storage.rows.is_empty()
    }
    pub fn mz(&self) -> &'a [f64] {
        &self.storage.mz
    }
    pub fn intensity(&self) -> &'a [f32] {
        &self.storage.intensity
    }
    pub fn get(&self, index: usize) -> Option<PeakRowView<'a>> {
        (index < self.len()).then_some(PeakRowView {
            peaks: *self,
            index,
        })
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = PeakRowView<'a>> + use<'a> {
        let peaks = *self;
        (0..peaks.len()).map(move |index| PeakRowView { peaks, index })
    }
}
#[derive(Debug)]
pub struct PeakIter<'a> {
    peaks: PeaksView<'a>,
    rows: Range<usize>,
}
impl<'a> Iterator for PeakIter<'a> {
    type Item = PeakRowView<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        self.rows.next().map(|index| PeakRowView {
            peaks: self.peaks,
            index,
        })
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.rows.size_hint()
    }
}
impl ExactSizeIterator for PeakIter<'_> {}
impl<'a> IntoIterator for PeaksView<'a> {
    type Item = PeakRowView<'a>;
    type IntoIter = PeakIter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        PeakIter {
            peaks: self,
            rows: 0..self.len(),
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct PeakRowView<'a> {
    peaks: PeaksView<'a>,
    index: usize,
}
impl<'a> PeakRowView<'a> {
    pub fn mz(&self) -> f64 {
        self.peaks.storage.mz[self.index]
    }
    pub fn intensity(&self) -> f32 {
        self.peaks.storage.intensity[self.index]
    }
    pub fn source_row(&self) -> usize {
        self.index
    }
    pub fn source_position(&self) -> SourcePosition {
        self.peaks.storage.rows[self.index].position
    }
    fn field(&self, range: &Option<Range<usize>>) -> SuppliedField<'a> {
        range.as_ref().map_or(SuppliedField::Absent, |r| {
            SuppliedField::Present(&self.peaks.record.raw.text[r.clone()])
        })
    }
    pub fn annotation_field(&self) -> SuppliedField<'a> {
        self.field(&self.peaks.storage.rows[self.index].annotation)
    }
    pub fn aggregation_field(&self) -> SuppliedField<'a> {
        self.field(&self.peaks.storage.rows[self.index].aggregation)
    }
    pub fn extra_columns(&self) -> impl Iterator<Item = &'a str> {
        let text = &self.peaks.record.raw.text;
        self.peaks.storage.rows[self.index]
            .extras
            .as_ref()
            .into_iter()
            .flat_map(move |r| text[r.clone()].split('\t'))
    }
}

#[derive(Debug, Default)]
pub(super) struct AnnotationStorage {
    alternatives: Vec<AnnotationOutcome>,
    rows: Vec<Range<usize>>,
    diagnostics: Vec<AnnotationDiagnostic>,
}
impl Clear for AnnotationStorage {
    fn clear_reuse(&mut self) {
        self.alternatives.clear();
        self.rows.clear();
        self.diagnostics.clear();
    }
}
#[derive(Debug)]
pub struct AnnotationOutcome {
    pub reported: PeakAnnotation,
    pub fragment: Option<Fragment<OutputMolecularFormula>>,
    pub diagnostic: Option<usize>,
}
#[derive(Debug)]
pub struct AnnotationDiagnostic {
    pub source_row: usize,
    pub alternative: Option<usize>,
    pub error: RecordError,
}
#[derive(Clone, Copy, Debug)]
pub struct AnnotationReportView<'a> {
    peaks: PeaksView<'a>,
    storage: &'a AnnotationStorage,
}
impl<'a> AnnotationReportView<'a> {
    pub fn len(&self) -> usize {
        self.storage.rows.len()
    }
    pub fn is_empty(&self) -> bool {
        self.storage.rows.is_empty()
    }
    pub fn get(&self, row: usize) -> Option<AnnotationRowView<'a>> {
        (row < self.len()).then_some(AnnotationRowView { report: *self, row })
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = AnnotationRowView<'a>> + use<'a> {
        let report = *self;
        (0..self.len()).map(move |row| AnnotationRowView { report, row })
    }
    pub fn diagnostics(&self) -> impl Iterator<Item = &'a AnnotationDiagnostic> + use<'a> {
        self.storage.diagnostics.iter()
    }
    pub fn has_diagnostics(&self) -> bool {
        !self.storage.diagnostics.is_empty()
    }
}
#[derive(Clone, Copy, Debug)]
pub struct AnnotationRowView<'a> {
    report: AnnotationReportView<'a>,
    row: usize,
}
impl<'a> AnnotationRowView<'a> {
    pub fn source_row(&self) -> usize {
        self.row
    }
    pub fn source_position(&self) -> SourcePosition {
        self.report.peaks.storage.rows[self.row].position
    }
    pub fn annotation_field(&self) -> SuppliedField<'a> {
        self.report.peaks.get(self.row).unwrap().annotation_field()
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &'a AnnotationOutcome> + use<'a> {
        self.report.storage.alternatives[self.report.storage.rows[self.row].clone()].iter()
    }
}

impl SpectrumRecord<'_> {
    pub fn peaks(&self) -> Result<PeaksView<'_>, RecordError> {
        self.require()?;
        let storage = self.peaks.get(|out| {
            let mut offset = self.raw.peaks.start;
            let first_line = self.raw.line_ends.partition_point(|end| *end <= offset);
            for (number, end) in self.raw.line_ends[first_line..].iter().copied().enumerate() {
                let full = &self.raw.text[offset..end];
                let line = full.trim_end_matches(['\r', '\n']);
                let pos = SourcePosition {
                    source: self.origin.source,
                    line: self.origin.line + first_line as u64 + number as u64,
                    byte_offset: self.origin.byte_offset + offset as u64,
                };
                if !line.trim().is_empty() && !line.starts_with('#') {
                    let mut fields = line.split('\t');
                    let mz_text = fields.next().unwrap_or_default();
                    let intensity_text = fields.next().ok_or_else(|| {
                        self.error(RecordErrorKind::Malformed, "missing peak intensity", pos)
                    })?;
                    let mz =
                        mz_text.trim().parse::<f64>().ok().filter(|n| n.is_finite()).ok_or_else(
                            || self.error(RecordErrorKind::Malformed, "invalid peak m/z", pos),
                        )?;
                    let intensity = intensity_text
                        .trim()
                        .parse::<f32>()
                        .ok()
                        .filter(|n| n.is_finite())
                        .ok_or_else(|| {
                            self.error(RecordErrorKind::Malformed, "invalid peak intensity", pos)
                        })?;
                    let mut at = offset + mz_text.len() + 1 + intensity_text.len() + 1;
                    let annotation = fields.next().map(|field| {
                        let r = at..at + field.len();
                        at += field.len() + 1;
                        r
                    });
                    let aggregation = fields.next().map(|field| {
                        let r = at..at + field.len();
                        at += field.len() + 1;
                        r
                    });
                    let extras = (at <= offset + line.len()).then_some(at..offset + line.len());
                    out.mz.push(mz);
                    out.intensity.push(intensity);
                    out.rows.push(Row {
                        annotation,
                        aggregation,
                        extras,
                        position: pos,
                    });
                }
                offset += full.len();
            }
            Ok(())
        })?;
        Ok(PeaksView {
            record: self,
            storage,
        })
    }

    pub fn analytes(&self) -> Result<&[DecodedAnalyte], RecordError> {
        self.require()?;
        self.analytes
            .get(|out| {
                for scope in self.analyte_scopes() {
                    let ScopeId::Analyte(id) = scope.id() else {
                        unreachable!()
                    };
                    if out.as_slice().iter().any(|a| a.id.get() == id) {
                        return Err(self.error(
                            RecordErrorKind::Conflict,
                            "duplicate analyte ID",
                            self.origin,
                        ));
                    }
                    let attrs = scope.attributes()?;
                    if out.active == out.slots.len() {
                        out.slots.push(DecodedAnalyte {
                            id: std::num::NonZeroU32::new(id).unwrap(),
                            target: AnalyteTarget::default(),
                            scratch: mzcore::sequence::ProFormaScratch::default(),
                            formulas: Reusable::default(),
                            spare_ion: PeptidoformIon::default(),
                            spare_formula: MolecularFormula::default(),
                        });
                    }
                    let analyte = &mut out.slots[out.active];
                    analyte.id = std::num::NonZeroU32::new(id).unwrap();
                    if let Some(a) = singleton(attrs, curie!(MS:1003270))? {
                        analyte.select_target(TargetKind::Peptidoform);
                        let AnalyteTarget::PeptidoformIon(ion) = &mut analyte.target else {
                            unreachable!()
                        };
                        analyte
                            .scratch
                            .parse_into(ion, a.raw_value(), self.context.ontologies)
                            .map_err(|e| {
                                self.error(
                                    RecordErrorKind::Malformed,
                                    format!("invalid ProForma: {e:?}"),
                                    a.origin().position,
                                )
                            })?;
                    } else if let Some(a) = singleton(attrs, curie!(MS:1000866))? {
                        analyte.select_target(TargetKind::Formula);
                        let AnalyteTarget::MolecularFormula(formula) = &mut analyte.target else {
                            unreachable!()
                        };
                        formula.pro_forma_into::<false, false>(a.raw_value()).map_err(|e| {
                            self.error(
                                RecordErrorKind::Malformed,
                                e.to_string(),
                                a.origin().position,
                            )
                        })?;
                    } else {
                        analyte.select_target(TargetKind::Unknown);
                    }
                    if let Some(charge) = singleton(attrs, curie!(MS:1000041))? {
                        let value = charge.raw_value().parse::<isize>().map_err(|_| {
                            self.error(
                                RecordErrorKind::Malformed,
                                "invalid analyte charge",
                                charge.origin().position,
                            )
                        })?;
                        let charge = Charge::new::<mzcore::system::e>(value);
                        if let AnalyteTarget::PeptidoformIon(ion) = &mut analyte.target {
                            analyte.scratch.set_charge(ion, charge);
                        } else {
                            analyte.target.set_charge(charge);
                        }
                    }
                    out.active += 1;
                }
                Ok(())
            })
            .map(DecodedAnalytes::as_slice)
    }

    pub fn interpretations(&self) -> Result<&[DecodedInterpretation], RecordError> {
        self.require()?;
        self.interpretations
            .get(|out| {
                for (scope_index, scope) in self.interpretation_scopes().enumerate() {
                    let ScopeId::Interpretation(id) = scope.id() else {
                        unreachable!()
                    };
                    let attrs = scope.attributes()?;
                    let probability =
                        singleton(attrs, curie!(MS:1002357))?.map(|a| a.to_f64()).transpose()?;
                    if out.active == out.slots.len() {
                        out.slots.push(DecodedInterpretation::default());
                    }
                    let interpretation = &mut out.slots[out.active];
                    interpretation.id = id;
                    interpretation.scope_index = scope_index;
                    interpretation.probability = probability;
                    interpretation.analyte_refs.clear();
                    for a in attrs.by_accession(curie!(MS:1003163)) {
                        for id in a.raw_value().split(',') {
                            interpretation.analyte_refs.push(id.trim().parse().map_err(|_| {
                                self.error(
                                    RecordErrorKind::Malformed,
                                    "invalid analyte reference",
                                    a.origin().position,
                                )
                            })?);
                        }
                    }
                    out.active += 1;
                }
                Ok(())
            })
            .map(DecodedInterpretations::as_slice)
    }

    pub fn resolved_annotations(&self) -> Result<AnnotationReportView<'_>, RecordError> {
        let peaks = self.peaks()?;
        let storage = self.annotations.get(|out| {
            let analytes = self.analytes(); // Soft dependency: never propagate with ?.
            for row in peaks.iter() {
                let start = out.alternatives.len();
                if let SuppliedField::Present(text) = row.annotation_field() {
                    if !text.trim().is_empty() && text.trim() != "?" {
                        match PeakAnnotation::parse_all(text.trim(), self.context.ontologies) {
                            Err(e) => out.diagnostics.push(AnnotationDiagnostic {
                                source_row: row.index,
                                alternative: None,
                                error: self.error(
                                    RecordErrorKind::Malformed,
                                    e.to_string(),
                                    row.source_position(),
                                ),
                            }),
                            Ok(reported) => {
                                for (alternative, annotation) in reported.into_iter().enumerate() {
                                    // The reported chemical expression and derived fragment are distinct
                                    // owned chemical objects. Raw text and analyte targets are never cloned.
                                    let resolved = match &analytes {
                                        Ok(targets) => annotation
                                            .clone()
                                            .into_fragment_with(
                                                |id| {
                                                    targets
                                                        .iter()
                                                        .find(|a| a.id.get() == id)
                                                        .map(|a| &a.target)
                                                },
                                                &Context::default().lines(0, text),
                                            )
                                            .map_err(|e| {
                                                self.error(
                                                    RecordErrorKind::Unsupported,
                                                    e.to_string(),
                                                    row.source_position(),
                                                )
                                            }),
                                        Err(e) => Err(self.error(
                                            RecordErrorKind::Unresolved,
                                            format!("annotation chemistry blocked: {e}"),
                                            row.source_position(),
                                        )),
                                    };
                                    let (fragment, diagnostic) = match resolved {
                                        Ok(f) => (Some(f), None),
                                        Err(error) => {
                                            let index = out.diagnostics.len();
                                            out.diagnostics.push(AnnotationDiagnostic {
                                                source_row: row.index,
                                                alternative: Some(alternative),
                                                error,
                                            });
                                            (None, Some(index))
                                        }
                                    };
                                    out.alternatives.push(AnnotationOutcome {
                                        reported: annotation,
                                        fragment,
                                        diagnostic,
                                    });
                                }
                            }
                        }
                    }
                }
                out.rows.push(start..out.alternatives.len());
            }
            Ok(())
        })?;
        Ok(AnnotationReportView { peaks, storage })
    }

    pub fn validate(&self) -> Result<(), RecordError> {
        self.validate_sections()?;
        self.project_description()?;
        for analyte in self.analytes()? {
            self.project_analyte(analyte.id, AnalyteTarget::default())?;
        }
        Ok(())
    }
    fn validate_sections(&self) -> Result<(), RecordError> {
        self.require()?;
        for scope in self.scopes() {
            for attribute in scope.attributes()?.iter() {
                attribute.value().map_err(Clone::clone)?;
            }
            if matches!(scope.id(), ScopeId::InterpretationMember(_)) {
                return Err(self.error(
                    RecordErrorKind::Unsupported,
                    "interpretation members are retained but not decoded",
                    self.origin,
                ));
            }
            if let Some(error) = scope.diagnostics()?.next() {
                return Err(error.clone());
            }
        }
        let analytes = self.analytes()?;
        for interpretation in self.interpretations()? {
            if interpretation
                .analyte_refs
                .iter()
                .any(|id| !analytes.iter().any(|a| a.id.get() == *id))
            {
                return Err(self.error(
                    RecordErrorKind::Unresolved,
                    "interpretation references missing analyte",
                    self.origin,
                ));
            }
        }
        if let Some(diagnostic) = self.resolved_annotations()?.diagnostics().next() {
            return Err(diagnostic.error.clone());
        }
        Ok(())
    }

    /// Explicit owned projection. Supplied evidence remains in this record.
    pub fn materialize(&self) -> Result<AnnotatedSpectrum<OutputMolecularFormula>, RecordError> {
        self.validate_sections()?;
        let mut spec = self.project_description()?;
        spec.analytes = self
            .analytes()?
            .iter()
            .map(|a| self.project_analyte(a.id, a.target.clone()))
            .collect::<Result<_, _>>()?;
        spec.interpretations = self
            .interpretations()?
            .iter()
            .map(|i| {
                self.project_interpretation(
                    i.scope_index,
                    i.id,
                    i.probability,
                    i.analyte_refs.clone(),
                )
            })
            .collect::<Result<_, _>>()?;
        let report = self.resolved_annotations()?;
        for peak in report.peaks.iter() {
            let annotations = report
                .get(peak.index)
                .unwrap()
                .iter()
                .filter_map(|a| a.fragment.clone())
                .collect();
            let aggregations = match peak.aggregation_field() {
                SuppliedField::Absent => Vec::new(),
                SuppliedField::Present(s) => s.split(',').map(str::to_owned).collect(),
            };
            spec.peaks.push(AnnotatedPeak::new(
                MassOverCharge::new::<mzcore::system::thomson>(peak.mz()),
                peak.intensity(),
                peak.index as u32,
                annotations,
                aggregations,
            ));
        }
        Ok(spec)
    }
    fn project_interpretation(
        &self,
        scope_index: usize,
        id: u32,
        probability: Option<f64>,
        analyte_refs: Vec<u32>,
    ) -> Result<Interpretation, RecordError> {
        let scope = self.interpretation_scopes().nth(scope_index).unwrap();
        Ok(Interpretation {
            id,
            probability,
            analyte_refs,
            attributes: flatten(owned_groups(scope.attributes()?)?),
            ..Interpretation::default()
        })
    }
    fn project_analyte(
        &self,
        id: std::num::NonZeroU32,
        target: AnalyteTarget,
    ) -> Result<Analyte, RecordError> {
        let scope = self
            .analyte_scopes()
            .find(|scope| scope.id() == ScopeId::Analyte(id.get()))
            .unwrap();
        let attrs = scope.attributes()?;
        let mut analyte = Analyte::new(id, target);
        let owned = owned_groups(attrs)?;
        for group in owned.values() {
            let mut protein = ProteinDescription::default();
            for (a, context) in group {
                if !protein.populate_from_attribute(a, context).map_err(|e| {
                    self.error(RecordErrorKind::Malformed, e.to_string(), self.origin)
                })? {
                    if ![curie!(MS:1003270), curie!(MS:1000866), curie!(MS:1000041)]
                        .contains(&a.name.accession)
                    {
                        analyte.params.push(a.clone().into());
                    }
                }
            }
            if !protein.is_empty() {
                analyte.proteins.push(protein);
            }
        }
        Ok(analyte)
    }
    fn project_description(
        &self,
    ) -> Result<AnnotatedSpectrum<OutputMolecularFormula>, RecordError> {
        let mut spec = AnnotatedSpectrum::<OutputMolecularFormula> {
            key: self.key().unwrap(),
            attributes: vec![Vec::new()],
            ..Default::default()
        };
        spec.description.ms_level = 2;
        spec.description.signal_continuity = mzdata::spectrum::SignalContinuity::Centroid;
        spec.description.acquisition.scans.push(ScanEvent::default());
        spec.description.acquisition.scans[0]
            .scan_windows
            .push(ScanWindow::default());
        spec.description.precursor.push(Precursor::default());
        spec.description.precursor[0].ions.push(SelectedIon::default());
        let groups = owned_groups(self.attributes()?)?;
        populate_spectrum_description_from_attributes(
            groups.iter(),
            &mut spec.description,
            &mut spec.attributes,
        )
        .map_err(|e| self.error(RecordErrorKind::Malformed, e.to_string(), self.origin))?;
        spec.description.precursor[0].activation._extract_methods_from_params();
        Ok(spec)
    }
    /// Consume the record, moving decoded analytes, interpretations and fragments.
    pub fn into_annotated(
        mut self,
    ) -> Result<AnnotatedSpectrum<OutputMolecularFormula>, RecordError> {
        self.validate_sections()?;
        let mut spec = self.project_description()?;
        spec.analytes = self
            .analytes
            .completed
            .take()
            .unwrap()
            .storage
            .into_active()
            .map(|a| self.project_analyte(a.id, a.target))
            .collect::<Result<_, _>>()?;
        spec.interpretations = self
            .interpretations
            .completed
            .take()
            .unwrap()
            .storage
            .into_active()
            .map(|i| {
                self.project_interpretation(i.scope_index, i.id, i.probability, i.analyte_refs)
            })
            .collect::<Result<_, _>>()?;
        let peaks = self.peaks.completed.take().unwrap().storage;
        let report = self.annotations.completed.take().unwrap().storage;
        let mut alternatives = report.alternatives.into_iter();
        for (index, (row, range)) in peaks.rows.into_iter().zip(report.rows).enumerate() {
            let annotations = alternatives
                .by_ref()
                .take(range.len())
                .filter_map(|a| a.fragment)
                .collect();
            let aggregations = row.aggregation.map_or_else(Vec::new, |r| {
                self.raw.text[r].split(',').map(str::to_owned).collect()
            });
            spec.peaks.push(AnnotatedPeak::new(
                MassOverCharge::new::<mzcore::system::thomson>(peaks.mz[index]),
                peaks.intensity[index],
                index as u32,
                annotations,
                aggregations,
            ));
        }
        Ok(spec)
    }
}

fn singleton(
    attrs: AttributeView<'_>,
    accession: Curie,
) -> Result<Option<AttributeOccurrence<'_>>, RecordError> {
    let mut values = attrs.by_accession(accession);
    let value = values.next();
    if let Some(other) = values.next() {
        return Err(RecordError::new(
            RecordErrorKind::Conflict,
            "multiple singleton values",
            other.origin().position,
        ));
    }
    Ok(value)
}

fn owned_groups(attrs: AttributeView<'_>) -> Result<super::super::AttributeGroups, RecordError> {
    let mut groups = super::super::AttributeGroups::new();
    groups.insert(None, Vec::new());
    for a in attrs.iter() {
        let line = format!("{}={}", a.name(), a.raw_value());
        let (_, attr, _) = Attribute::parse(&line).map_err(|e| {
            RecordError::new(
                RecordErrorKind::Malformed,
                e.to_string(),
                a.origin().position,
            )
        })?;
        let context = Context::default()
            .line_index(a.origin().position.line as u32)
            .lines(0, line)
            .to_owned();
        let group = a
            .group_id()
            .map(|g| {
                u32::try_from(g.number()).map_err(|_| {
                    RecordError::new(
                        RecordErrorKind::Malformed,
                        "group ID exceeds materialized representation",
                        a.origin().position,
                    )
                })
            })
            .transpose()?;
        groups.entry(group).or_default().push((attr, context));
    }
    Ok(groups)
}
fn flatten(groups: super::super::AttributeGroups) -> Attributes {
    groups
        .into_values()
        .map(|g| g.into_iter().map(|(a, _)| a).collect())
        .collect()
}
