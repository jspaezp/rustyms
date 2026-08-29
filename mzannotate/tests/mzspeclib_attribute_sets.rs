//! Test that library entries inherit the attribute sets they claim, per mzSpecLib v1.0 §4.1.11.
#![allow(clippy::missing_panics_doc)]
use std::collections::HashMap;

use context_error as _;
use indexmap as _;
use itertools as _;
use mzannotate::mzspeclib::MzSpecLibTextParser;
use mzcv as _;
use mzdata as _;
use num_traits as _;
use ordered_float as _;
use serde as _;
use serde_json as _;
use thin_vec as _;
use uom as _;

/// The five decoys in this library claim `<AttributeSet Spectrum=DECOY>`, which overrides the
/// `observed spectrum` origin type that the `all` set gives every entry.
const DECOYS: &[&str] = &[
    "AVC[Carbamidomethyl]ASFSLTHR/3 (HCD)",
    "VLTFPVANK/2 (HCD)",
    "APFLLGEITFTGHFFPK/3 (HCD)",
    "NLASQVGWGM[Oxidation]VHAR/2 (HCD)",
    "HFDPDGPGSTAVNNTR/2 (HCD)",
];

/// Attributes the mzSpecLib metadata rules mark as `Repeatable | False`, so a resolved spectrum
/// must carry at most one value for each.
const NOT_REPEATABLE: &[&str] = &[
    "MS:1003072", // spectrum origin type
    "MS:1000044", // dissociation method
    "MS:1000031", // instrument model
];

#[test]
fn claimed_attribute_sets_are_applied_and_override() {
    let spectra = MzSpecLibTextParser::open_file(
        "../data/fetal_brain_tiny_consensus_td.mzSpecLib.txt".as_ref(),
        &mzcore::ontology::STATIC_ONTOLOGIES,
    )
    .unwrap()
    .collect::<Result<Vec<_>, _>>()
    .unwrap();

    assert_eq!(spectra.len(), 10);

    for spectrum in &spectra {
        let origin: Vec<_> = spectrum
            .description
            .params
            .iter()
            .filter(|p| p.curie().is_some_and(|c| c.to_string() == "MS:1003072"))
            .map(|p| p.value.to_string())
            .collect();
        let expected = if DECOYS.contains(&spectrum.description.id.as_str()) {
            "MS:1003195|shuffle-and-reposition decoy spectrum"
        } else {
            "MS:1003073|observed spectrum"
        };
        assert_eq!(
            origin,
            vec![expected.to_string()],
            "wrong spectrum origin type for {}",
            spectrum.description.id
        );
    }
}

#[test]
fn non_repeatable_attributes_resolve_to_one_value() {
    for entry in std::fs::read_dir("../data").unwrap().flatten() {
        if !entry
            .path()
            .to_string_lossy()
            .to_ascii_lowercase()
            .ends_with(".mzspeclib.txt")
        {
            continue;
        }

        for spectrum in MzSpecLibTextParser::open_file(
            &entry.path(),
            &mzcore::ontology::STATIC_ONTOLOGIES,
        )
        .unwrap()
        .flatten()
        {
            let mut counts: HashMap<String, usize> = HashMap::new();
            for param in &spectrum.description.params {
                if let Some(curie) = param.curie() {
                    *counts.entry(curie.to_string()).or_default() += 1;
                }
            }
            for accession in NOT_REPEATABLE {
                assert!(
                    counts.get(*accession).copied().unwrap_or(0) <= 1,
                    "{} carries {} values for {accession} in {}",
                    spectrum.description.id,
                    counts[*accession],
                    entry.path().display(),
                );
            }
        }
    }
}
