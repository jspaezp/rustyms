//! Reusable storage for common linear ProForma ions. Complex syntax delegates to
//! the full parser, preserving its results and diagnostics.
use crate::{
    ParserResult,
    chemistry::{MolecularCharge, MolecularFormula},
    ontology::{Ontologies, Ontology},
    sequence::{
        AminoAcid, CheckedAminoAcid, Linked, Modification, Peptidoform, PeptidoformIon,
        PlacementRule, SequenceElement, SequencePosition, SimpleModification,
        SimpleModificationInner,
    },
    system::isize::Charge,
};
use context_error::BasicKind;
use thin_vec::ThinVec;

/// Scratch storage for repeatedly decoding ProForma into the same ion.
/// Linear uppercase sequences with numeric UNIMOD/MOD modifications and optional
/// integer proton charge reuse residue/modification/charge storage. Other syntax
/// retains full parser semantics through an allocating fallback.
#[derive(Debug, Default)]
pub struct ProFormaScratch {
    residues: Vec<SequenceElement<Linked>>,
    modifications: Vec<(usize, SimpleModification)>,
    charge: Option<MolecularCharge>,
    spare: PeptidoformIon,
    simple_output: bool,
}
impl ProFormaScratch {
    /// Replace `output` with a parsed ion. On failure its contents are unspecified;
    /// the scratch allocations remain reusable. Returns ordinary parser warnings.
    pub fn parse_into<'a>(
        &mut self,
        output: &mut PeptidoformIon,
        text: &'a str,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, (), BasicKind> {
        if let Some((sequence_end, length, charge)) = self.scan(text, ontologies) {
            if !self.simple_output {
                std::mem::swap(output, &mut self.spare);
                self.simple_output = true;
            }
            if output.peptidoforms.is_empty() {
                output.peptidoforms.push(Peptidoform::default());
            }
            output.peptidoforms.truncate(1);
            output.name_mut().clear();
            let peptide = &mut output.peptidoforms[0];
            peptide.prepare_linear_reuse();
            let sequence = peptide.sequence_mut();
            while sequence.len() > length {
                self.residues.push(sequence.pop().unwrap());
            }
            while sequence.len() < length {
                sequence.push(self.residues.pop().unwrap_or_default());
            }
            let mut index = 0;
            let mut depth = false;
            for byte in text.as_bytes()[..sequence_end].iter().copied() {
                match byte {
                    b'[' => depth = true,
                    b']' => depth = false,
                    _ if !depth => {
                        let residue = &mut sequence[index];
                        residue.aminoacid =
                            CheckedAminoAcid::new(AminoAcid::try_from(byte).unwrap());
                        residue.modifications.clear();
                        residue.ambiguous = None;
                        index += 1;
                    }
                    _ => (),
                }
            }
            for (index, modification) in &self.modifications {
                sequence[*index]
                    .modifications
                    .push(Modification::Simple(modification.clone()));
            }
            self.replace_charge(output, charge);
            return Ok(((), Vec::new()));
        }
        if self.simple_output {
            std::mem::swap(output, &mut self.spare);
            self.simple_output = false;
        }
        let (parsed, warnings) = PeptidoformIon::pro_forma(text, ontologies)?;
        *output = parsed;
        Ok(((), warnings))
    }

    /// Apply a reported charge, preserving existing carrier chemistry if its
    /// charge already matches (the same rule as mzSpecLib analyte projection).
    pub fn set_charge(&mut self, output: &mut PeptidoformIon, charge: Charge) {
        if output.get_charge_carriers().is_some_and(|c| c.charge() == charge) {
            return;
        }
        if self.simple_output {
            self.replace_charge(output, Some(charge.value));
        } else {
            output.set_charge_carriers(Some(MolecularCharge::proton(charge)));
        }
    }

    fn replace_charge(&mut self, output: &mut PeptidoformIon, charge: Option<isize>) {
        let peptide = &mut output.peptidoforms[0];
        let mut carriers = peptide
            .take_charge_for_reuse()
            .filter(|c| {
                c.charge_carriers.len() == 1 && {
                    let formula = &c.charge_carriers[0].1;
                    formula.additional_mass() == 0.0
                        && formula.elements()
                            == &[
                                (crate::chemistry::Element::Electron, None, -1),
                                (crate::chemistry::Element::H, None, 1),
                            ]
                }
            })
            .or_else(|| self.charge.take());
        match charge {
            Some(0) => {
                self.charge = carriers;
                peptide.restore_charge_for_reuse(Some(MolecularCharge {
                    charge_carriers: ThinVec::new(),
                }));
            }
            Some(charge) => {
                let carrier = carriers.get_or_insert_with(|| {
                    MolecularCharge::proton(Charge::new::<crate::system::e>(1))
                });
                carrier.charge_carriers[0].0 = charge;
                peptide.restore_charge_for_reuse(carriers);
            }
            None => self.charge = carriers,
        }
    }

    fn scan(
        &mut self,
        text: &str,
        ontologies: &Ontologies,
    ) -> Option<(usize, usize, Option<isize>)> {
        self.modifications.clear();
        let (sequence, charge) = match text.split_once('/') {
            Some((sequence, charge)) => {
                let charge = charge.parse::<isize>().ok()?;
                if charge == 0 {
                    return None;
                }
                (sequence, Some(charge))
            }
            None => (text, None),
        };
        let mut at = 0;
        let mut count = 0;
        while at < sequence.len() {
            let byte = sequence.as_bytes()[at];
            if !byte.is_ascii_uppercase() {
                return None;
            }
            AminoAcid::try_from(byte).ok()?;
            at += 1;
            while sequence.as_bytes().get(at) == Some(&b'[') {
                let end = sequence[at + 1..].find(']')? + at + 1;
                let (prefix, id) = sequence[at + 1..end].split_once(':')?;
                let ontology = match prefix {
                    "UNIMOD" => Ontology::Unimod,
                    "MOD" => Ontology::Psimod,
                    _ => return None,
                };
                let id = id.parse::<u32>().ok()?;
                let modification = ontologies.get_by_index(ontology, &id.into())?;
                let SimpleModificationInner::Database { .. } = modification.as_ref() else {
                    return None;
                };
                self.modifications.push((count, modification));
                at = end + 1;
            }
            count += 1;
        }
        if count == 0 {
            return None;
        }
        // Use the same placement predicates as the full parser. Nonmatching
        // placements take its diagnostic path, including non-fatal warnings.
        let mut index = 0;
        let mut depth = false;
        for byte in sequence.bytes() {
            match byte {
                b'[' => depth = true,
                b']' => depth = false,
                _ if !depth => {
                    let residue = SequenceElement::<Linked>::new(
                        CheckedAminoAcid::new(AminoAcid::try_from(byte).ok()?),
                        None,
                    );
                    for (_, modification) in self.modifications.iter().filter(|(i, _)| *i == index)
                    {
                        if let SimpleModificationInner::Database { specificities, .. } =
                            modification.as_ref()
                        {
                            if !specificities.is_empty()
                                && !specificities.iter().any(|(rules, ..)| {
                                    PlacementRule::any_possible(
                                        rules,
                                        &residue,
                                        SequencePosition::Index(index, count),
                                    )
                                })
                            {
                                return None;
                            }
                        }
                    }
                    index += 1;
                }
                _ => (),
            }
        }
        Some((sequence.len(), count, charge))
    }
}

/// Reusable formula outputs for a peptidoform ion. Common unambiguous linear
/// formulas are accumulated directly; general chemistry uses the established
/// resolver and copies its outputs into retained slots.
#[derive(Debug, Default)]
pub struct FormulaBuffer {
    values: Vec<MolecularFormula>,
    active: usize,
}
impl FormulaBuffer {
    /// Invalidate previous results without releasing their storage.
    pub fn clear(&mut self) {
        self.active = 0;
    }
    /// Borrow the current results.
    pub fn as_slice(&self) -> &[MolecularFormula] {
        &self.values[..self.active]
    }
    /// Copy a reported formula into retained storage.
    pub fn set_formula(&mut self, formula: &MolecularFormula) {
        if self.values.is_empty() {
            self.values.push(MolecularFormula::default());
        }
        self.values[0].clear();
        self.values[0] += formula;
        self.active = 1;
    }
    /// Calculate formulas, reusing output capacity across calls.
    pub fn calculate(&mut self, ion: &PeptidoformIon) -> &[MolecularFormula] {
        use crate::chemistry::AmbiguousMolecule;
        if self.values.is_empty() {
            self.values.push(MolecularFormula::default());
        }
        self.active = 0;
        if ion.peptidoforms.len() == 1
            && ion.peptidoforms[0].formula_into_reuse(&mut self.values[0])
        {
            self.active = 1;
        } else {
            for formula in ion.formulas().iter() {
                if self.active == self.values.len() {
                    self.values.push(MolecularFormula::default());
                }
                self.values[self.active].clear();
                self.values[self.active] += formula;
                self.active += 1;
            }
        }
        self.as_slice()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        chemistry::{AmbiguousMolecule, MassOutputType},
        ontology::STATIC_ONTOLOGIES,
    };

    #[test]
    fn reuse_matches_full_parser_across_shapes_failures_and_fallbacks() {
        let mut scratch = ProFormaScratch::default();
        let mut output = PeptidoformIon::default();
        let mut formulas = FormulaBuffer::default();
        for text in [
            "AC[UNIMOD:4]M[UNIMOD:35]/2",
            "A",
            "PEPTIDE/3",
            "M[UNIMOD:35]M[UNIMOD:35]/-2",
            "[UNIMOD:1]-PEPTIDE/2",
            "<13C>PEPTIDE/2",
            "BZX/2",
            "PEP[+15.5]TIDE/2",
            "A[UNIMOD:4]/2",
            "AC[UNIMOD:4]M[UNIMOD:35]/2",
            "???",
            "M[UNIMOD:35]/2",
            "C[X:Disulfide#xl1]PEPC[#xl1]/2",
            "PEPTIDE/2[+2Na+]",
            "PEPTIDE",
            "A/0",
            "A[UNIMOD:999999999]",
            "A[UNIMOD:4",
            "AC[UNIMOD:4]/2",
            "M[MOD:00719]/2",
        ]
        .into_iter()
        .cycle()
        .take(60)
        {
            let expected = PeptidoformIon::pro_forma(text, &STATIC_ONTOLOGIES);
            let observed = scratch.parse_into(&mut output, text, &STATIC_ONTOLOGIES);
            match (expected, observed) {
                (Ok((expected, warnings)), Ok(((), observed_warnings))) => {
                    assert_eq!(output, expected, "{text}");
                    assert_eq!(
                        format!("{warnings:?}"),
                        format!("{observed_warnings:?}"),
                        "{text}"
                    );
                    let expected = expected.formulas();
                    let observed = formulas.calculate(&output);
                    assert_eq!(observed, &*expected, "{text}");
                    for (a, b) in observed.iter().zip(expected.iter()) {
                        assert_eq!(a.labels(), b.labels(), "{text}");
                    }
                }
                (Err(expected), Err(observed)) => {
                    assert_eq!(format!("{expected:?}"), format!("{observed:?}"), "{text}")
                }
                (expected, observed) => panic!("{text}: {expected:?} != {observed:?}"),
            }
        }
    }

    #[test]
    fn formula_labels_keep_general_resolver_order() {
        use crate::chemistry::AmbiguousLabel;
        let mut ion = PeptidoformIon::pro_forma("PEPTIDE", &STATIC_ONTOLOGIES).unwrap().0;
        let labelled = |index| Modification::Simple(std::sync::Arc::new(SimpleModificationInner::Formula(
            crate::molecular_formula!(C 1).with_label(AmbiguousLabel::AminoAcid {
                option: AminoAcid::Alanine, sequence_index: index, peptidoform_index: 0, peptidoform_ion_index: 0,
            })
        )));
        ion.peptidoforms_mut()[0].set_n_term(vec![labelled(0)]);
        ion.peptidoforms_mut()[0].set_c_term(vec![labelled(6)]);
        let expected = ion.formulas();
        let mut buffer = FormulaBuffer::default();
        let observed = buffer.calculate(&ion);
        assert_eq!(observed, &*expected);
        for (a, b) in observed.iter().zip(expected.iter()) { assert_eq!(a.labels(), b.labels()); }
    }

    #[test]
    fn changing_charge_and_external_output_does_not_keep_stale_chemistry() {
        let mut scratch = ProFormaScratch::default();
        let mut output = PeptidoformIon::default();
        scratch.parse_into(&mut output, "PEPTIDE/2", &STATIC_ONTOLOGIES).unwrap();
        for charge in [0, 3, -2, 0, 2] {
            scratch.set_charge(&mut output, Charge::new::<crate::system::e>(charge));
            assert_eq!(
                output.get_charge_carriers().unwrap(),
                &MolecularCharge::proton(Charge::new::<crate::system::e>(charge))
            );
        }
        output = PeptidoformIon::pro_forma("<13C>[UNIMOD:1]-PEPTIDE/2", &STATIC_ONTOLOGIES)
            .unwrap()
            .0;
        scratch
            .parse_into(&mut output, "AC[UNIMOD:4]/3", &STATIC_ONTOLOGIES)
            .unwrap();
        assert_eq!(
            output,
            PeptidoformIon::pro_forma("AC[UNIMOD:4]/3", &STATIC_ONTOLOGIES).unwrap().0
        );
    }
}
