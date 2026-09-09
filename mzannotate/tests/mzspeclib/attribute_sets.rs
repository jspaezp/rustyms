//! Parser regressions for attribute-set inheritance, precedence, and group scope.
//!
//! mzSpecLib v1.0 specification (May 23, 2025), section 4.1.4, pp. 13-14,
//! "Attribute set definitions": `all` applies implicitly; named sets apply via
//! MS:1003212; later claims override earlier ones; entry values override sets.
//! Repeated instances of a term are inherited and replaced together.
//! https://github.com/HUPO-PSI/mzSpecLib/blob/506791706f41e95e0c2a1d4de15e49f25b59c760/specification/mzSpecLib_specification_v1.0.pdf#page=14
//!
//! Section 4.1.11, p. 21, specifies replacement rather than accumulation.
//! Section 4.1.12 supplies worked examples: claim order (Example 1, pp. 21-23),
//! repeated terms (Example 2, p. 23), and grouped claims (Example 3, pp. 23-24).
use mzannotate::{
    mzspeclib::{Attribute, Attributes, MzSpecLibTextParser},
    spectrum::AnnotatedSpectrum,
};
use mzcore::chemistry::OutputMolecularFormula;
use mzcv::curie;

fn parse_spectrum(text: &str) -> AnnotatedSpectrum<OutputMolecularFormula> {
    let mut parser =
        MzSpecLibTextParser::open(text.as_bytes(), None, &mzcore::ontology::STATIC_ONTOLOGIES)
            .expect("valid library header");
    let spectrum = parser
        .next()
        .expect("fixture contains one spectrum")
        .expect("valid spectrum");
    assert!(
        parser.next().is_none(),
        "fixture contains exactly one spectrum"
    );
    spectrum
}

fn action_values(attributes: &[Attribute]) -> Vec<String> {
    attributes
        .iter()
        .filter(|attribute| attribute.name.accession == curie!(MS:1000543))
        .map(|attribute| attribute.value.to_string())
        .collect()
}

fn interpretation_attributes(claims_and_attributes: &str) -> Attributes {
    // Synthetic interpretation fixture exercising the type-independent rules in
    // sections 4.1.4 and 4.1.12; not copied from an upstream interpretation example.
    // Declaration order deliberately differs from the claim order tested below.
    let text = format!(
        "<mzSpecLib>
MS:1003186|library format version=1.0
<AttributeSet Interpretation=all>
MS:1000543|data processing action=default first
MS:1000543|data processing action=default second
<AttributeSet Interpretation=later>
MS:1000543|data processing action=later named
<AttributeSet Interpretation=named>
MS:1000543|data processing action=named first
MS:1000543|data processing action=named second
<Spectrum=1>
<Interpretation=1>
{claims_and_attributes}
<Peaks>"
    );
    let spectrum = parse_spectrum(&text);
    assert_eq!(spectrum.interpretations.len(), 1);
    spectrum.interpretations[0].attributes.clone()
}

#[test]
fn spectrum_named_set_replaces_default_origin_in_description() {
    // Before: the named claim was missed, leaving only the default
    // MS:1003073|observed spectrum in description.params.
    // After: the claimed set overrides the default, leaving exactly one origin:
    // MS:1003195|shuffle-and-reposition decoy spectrum.
    // Reduced from the upstream SpectraST example: `all` supplies observed origin
    // (line 29), DECOY supplies decoy origin (lines 35-36), and line 246 claims it.
    // https://github.com/HUPO-PSI/mzSpecLib/blob/506791706f41e95e0c2a1d4de15e49f25b59c760/examples/SpectraST/fetal_brain_tiny_consensus_td.mzSpecLib.txt#L18-L36
    // https://github.com/HUPO-PSI/mzSpecLib/blob/506791706f41e95e0c2a1d4de15e49f25b59c760/examples/SpectraST/fetal_brain_tiny_consensus_td.mzSpecLib.txt#L246
    let spectrum = parse_spectrum(
        "<mzSpecLib>
MS:1003186|library format version=1.0
<AttributeSet Spectrum=all>
MS:1003072|spectrum origin type=MS:1003073|observed spectrum
<AttributeSet Spectrum=decoy>
MS:1003072|spectrum origin type=MS:1003195|shuffle-and-reposition decoy spectrum
<Spectrum=1>
MS:1003212|library attribute set name=decoy
<Peaks>",
    );
    let origins: Vec<_> = spectrum
        .description
        .params
        .iter()
        .filter(|param| param.curie().is_some_and(|curie| curie.to_string() == "MS:1003072"))
        .map(|param| param.value.to_string())
        .collect();
    assert_eq!(
        origins,
        ["MS:1003195|shuffle-and-reposition decoy spectrum"],
        "the claimed decoy origin must replace the default in parsed metadata"
    );
}

#[test]
fn interpretation_inherits_all_without_any_named_claim() {
    // Section 4.1.4, pp. 13-14: the reserved `all` set applies to every entry
    // of its type, without an explicit MS:1003212 claim.
    let attributes = interpretation_attributes("");
    assert_eq!(attributes.len(), 1, "only the ungrouped context is present");
    assert_eq!(
        action_values(&attributes[0]),
        ["default first", "default second"],
        "all applies without a claim and preserves repeated values"
    );
}

#[test]
fn named_set_replaces_every_default_value_and_preserves_its_repeats() {
    // Section 4.1.11, p. 21: named instances supersede the instances in `all`
    // together; multiple values within the winning context remain together.
    let attributes = interpretation_attributes("MS:1003212|library attribute set name=named");
    assert_eq!(
        action_values(&attributes[0]),
        ["named first", "named second"],
        "named values replace both defaults, rather than appending or deduplicating"
    );
}

#[test]
fn last_claimed_set_wins_regardless_of_header_declaration_order() {
    // Section 4.1.12, Example 1, pp. 22-23 explicitly distinguishes claim order
    // from the order in which sets are defined in the library header.
    let attributes = interpretation_attributes(
        "MS:1003212|library attribute set name=named
MS:1003212|library attribute set name=later",
    );
    assert_eq!(
        action_values(&attributes[0]),
        ["later named"],
        "the last claim replaces both earlier named values and all defaults"
    );
}

#[test]
fn entry_values_replace_every_inherited_value_and_preserve_its_repeats() {
    // Section 4.1.11, p. 21, and section 4.1.12, Example 2, p. 23:
    // local instances replace every inherited instance of the same term.
    let attributes = interpretation_attributes(
        "MS:1003212|library attribute set name=named
MS:1000543|data processing action=entry first
MS:1000543|data processing action=entry second",
    );
    assert_eq!(
        action_values(&attributes[0]),
        ["entry first", "entry second"],
        "both entry values survive, with no named or default values remaining"
    );
}

#[test]
fn grouped_claim_inherits_exact_value_only_inside_its_group() {
    // Section 4.1.12, Example 3, pp. 23-24: claiming a set inside a group
    // places its inherited attributes in that group.
    // Adapted to Interpretation from the upstream ANALYTE example, whose
    // HUMAN_TRYPTIC definition is at lines 30-35 and grouped claim at line 74.
    // https://github.com/HUPO-PSI/mzSpecLib/blob/506791706f41e95e0c2a1d4de15e49f25b59c760/examples/SpectraST/fetal_brain_tiny.mzSpecLib.txt#L30-L35
    // https://github.com/HUPO-PSI/mzSpecLib/blob/506791706f41e95e0c2a1d4de15e49f25b59c760/examples/SpectraST/fetal_brain_tiny.mzSpecLib.txt#L74-L77
    let spectrum = parse_spectrum(
        "<mzSpecLib>
MS:1003186|library format version=1.0
<AttributeSet Interpretation=human>
MS:1001469|taxonomy: scientific name=Homo sapiens
<Spectrum=1>
<Interpretation=1>
[2]MS:1003212|library attribute set name=human
<Peaks>",
    );
    assert_eq!(spectrum.interpretations.len(), 1);
    let attributes = &spectrum.interpretations[0].attributes;
    assert_eq!(
        attributes.len(),
        2,
        "ungrouped context plus the claiming group"
    );
    assert!(
        attributes[0].is_empty(),
        "grouped inheritance must not leak outside the group"
    );
    let taxonomy: Vec<_> = attributes[1]
        .iter()
        .filter(|attribute| attribute.name.accession == curie!(MS:1001469))
        .map(|attribute| attribute.value.to_string())
        .collect();
    assert_eq!(
        taxonomy,
        ["Homo sapiens"],
        "the claiming group inherits the exact taxonomy value once"
    );
}
