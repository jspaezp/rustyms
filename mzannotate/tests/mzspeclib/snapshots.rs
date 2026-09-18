//! Opt-in characterization capture; driven by scripts/mzspeclib_snapshots.py.
use std::{collections::BTreeMap, io::Cursor, path::Path};

use context_error::{BoxedError, ErrorKind, FullErrorContent};
use mzannotate::{
    mzspeclib::{
        AnalyteTarget, Attribute, AttributeValue, Attributes, CleaveAgent, LibraryHeader,
        MzSpecLibErrorKind, MzSpecLibTextParser, ProteinDescription,
    },
    spectrum::AnnotatedSpectrum,
};
use mzcore::{chemistry::OutputMolecularFormula, ontology::STATIC_ONTOLOGIES};
use serde_json::{Value, json};

type Spectrum = AnnotatedSpectrum<OutputMolecularFormula>;
type ParseResult = Result<Spectrum, BoxedError<'static, MzSpecLibErrorKind>>;

fn attribute(attribute: &Attribute) -> Value {
    let value = match &attribute.value {
        AttributeValue::Scalar(value) => json!({"scalar": value}),
        AttributeValue::List(values) => json!({"list": values}),
        AttributeValue::Term(term) => json!({"term": term}),
    };
    json!({"name": attribute.name, "value": value})
}

fn attributes(groups: &Attributes) -> Value {
    json!(
        groups
            .iter()
            .map(|group| group.iter().map(attribute).collect::<Vec<_>>())
            .collect::<Vec<_>>()
    )
}

fn header(header: &LibraryHeader) -> Value {
    let classes: BTreeMap<_, _> = header
        .attribute_classes
        .iter()
        .map(|(kind, sets)| {
            (
                kind.to_string(),
                sets.iter()
                    .map(|set| {
                        json!({
                            "id": set.id, "namespace": set.namespace.to_string(),
                            "attributes": attributes(&set.attributes),
                        })
                    })
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    json!({"version": header.format_version, "attributes": attributes(&header.attributes), "classes": classes})
}

fn protein(protein: &ProteinDescription) -> Value {
    let cleavage = match &protein.cleavage_agent {
        CleaveAgent::Unknown => json!({"unknown": true}),
        CleaveAgent::Name(name) => json!({"name": name}),
        CleaveAgent::Term(term) => json!({"term": term}),
    };
    json!({
        "accession": protein.accession, "name": protein.name,
        "database_name": protein.database_name, "database_version": protein.database_version,
        "cleavage_agent": cleavage, "description": protein.description,
        "flanking_sequences": protein.flanking_sequences,
        "species_accession": protein.species_accession,
        "species_common_name": protein.species_common_name,
        "species_scientific_name": protein.species_scientific_name,
        "missed_cleavages": protein.missed_cleavages, "enzymatic_termini": protein.enzymatic_termini,
        "set_names": protein.set_names,
        "attributes": protein.attributes.iter().map(attribute).collect::<Vec<_>>(),
    })
}

fn spectrum(spectrum: &Spectrum) -> Value {
    let analytes: Vec<_> = spectrum
        .analytes
        .iter()
        .map(|analyte| {
            let target = match &analyte.target {
                AnalyteTarget::Unknown(charge) => json!({"unknown": charge}),
                AnalyteTarget::MolecularFormula(formula) => json!({"formula": formula}),
                AnalyteTarget::PeptidoformIon(ion) => json!({"peptidoform_ion": ion}),
            };
            // read_analyte traverses HashMap groups to populate these collections.
            // Their relative order is already nondeterministic; retain multiplicity.
            let mut params: Vec<_> = analyte.params.iter().map(|p| json!(p)).collect();
            let mut proteins: Vec<_> = analyte.proteins.iter().map(protein).collect();
            params.sort_by_cached_key(Value::to_string);
            proteins.sort_by_cached_key(Value::to_string);
            json!({"id": analyte.id, "target": target, "params": params,
            "proteins": proteins})
        })
        .collect();
    let interpretations: Vec<_> = spectrum
        .interpretations
        .iter()
        .map(|interpretation| {
            let members: BTreeMap<_, _> = interpretation
                .members
                .iter()
                .map(|(id, attrs)| (*id, attributes(attrs)))
                .collect();
            json!({"id": interpretation.id, "probability": interpretation.probability,
            "attributes": attributes(&interpretation.attributes),
            "analyte_refs": interpretation.analyte_refs, "members": members})
        })
        .collect();
    json!({"key": spectrum.key, "description": spectrum.description,
        "attributes": attributes(&spectrum.attributes), "analytes": analytes,
        "interpretations": interpretations,
        "peaks": spectrum.peaks.iter().collect::<Vec<_>>()})
}

// Snapshot explicit error content, not a Debug dump or terminal rendering.
#[allow(single_use_lifetimes)]
fn error<'a, K: ErrorKind + std::fmt::Debug>(error: &impl FullErrorContent<'a, K>) -> Value {
    json!({"kind": format!("{:?}", error.get_kind()),
        "short": error.get_short_description(), "long": error.get_long_description(),
        "contexts": error.get_contexts(), "suggestions": error.get_suggestions(),
        "version": error.get_version(),
        "underlying": error.get_underlying_errors().iter().map(|e| self::error(e)).collect::<Vec<_>>()})
}

fn result(result: ParseResult) -> Value {
    match result {
        Ok(value) => json!({"ok": spectrum(&value)}),
        Err(value) => json!({"error": error(&value)}),
    }
}

// Indexed reads usually reproduce an already captured event. Refer to its full
// representation rather than storing the same peak arrays four more times.
fn lookup(value: Option<ParseResult>, events: &[Value]) -> Value {
    value.map_or(Value::Null, |value| {
        let value = result(value);
        events
            .iter()
            .position(|event| *event == value)
            .map_or(value, |index| json!({"sequential_event": index}))
    })
}

fn capture(root: &Path, relative: &str) -> Value {
    let bytes = std::fs::read(root.join(relative)).unwrap();
    let open = || {
        MzSpecLibTextParser::open(
            Cursor::new(bytes.as_slice()),
            Some(relative.into()),
            &STATIC_ONTOLOGIES,
        )
    };
    let mut parser = match open() {
        Ok(parser) => parser,
        Err(e) => return json!({"open_error": error(&e)}),
    };
    let header = header(parser.header());
    let mut events = Vec::new();
    let mut queries = Vec::new();
    // Fail the capture if a regression fails to make forward progress.
    let limit = bytes.iter().filter(|byte| **byte == b'\n').count() + 16;
    loop {
        let Some(event) = parser.next() else { break };
        assert!(events.len() < limit, "nonterminating reader: {relative}");
        if let Ok(spectrum) = &event {
            queries.push((spectrum.key, spectrum.description.id.clone()));
        }
        events.push(result(event));
    }
    let eof = [
        lookup(parser.next(), &events),
        lookup(parser.next(), &events),
    ];
    let read_next_after_eof = result(parser.read_next());

    let mut indexed = open().unwrap();
    let index = match indexed.build_index() {
        Err(e) => json!({"io_error": e.to_string()}),
        Ok(()) => {
            let len = indexed.len();
            let by_index: Vec<_> = (0..len)
                .map(|i| lookup(indexed.get_spectrum_by_index(i), &events))
                .collect();
            let by_key: Vec<_> = queries
                .iter()
                .map(|(key, _)| {
                    json!({"key": key,
                "result": lookup(indexed.get_spectrum_by_key(*key), &events)})
                })
                .collect();
            let by_name: Vec<_> = queries
                .iter()
                .map(|(_, name)| {
                    json!({"name": name,
                "result": lookup(indexed.get_spectrum_by_name(name), &events)})
                })
                .collect();
            // Exercise every literal scan-number declaration, including duplicates.
            let scans: Vec<usize> = String::from_utf8_lossy(&bytes)
                .lines()
                .filter_map(|line| {
                    line.strip_prefix("MS:1003057|scan number=")?.trim().parse().ok()
                })
                .collect();
            let by_scan: Vec<_> = scans
                .iter()
                .map(|scan| {
                    json!({"scan": scan,
                "result": lookup(indexed.get_spectrum_by_scan_number(*scan), &events)})
                })
                .collect();
            json!({"len": len, "is_empty": indexed.is_empty(), "by_index": by_index,
                "by_key": by_key, "by_name": by_name, "by_scan": by_scan,
                "missing_index": lookup(indexed.get_spectrum_by_index(len), &events),
                "next_after_lookups": lookup(indexed.next(), &events)})
        }
    };
    json!({"schema": 1, "fixture": relative, "header": header,
        "events": events, "eof": eof, "read_next_after_eof": read_next_after_eof, "indexed": index})
}

#[test]
#[ignore = "run through scripts/mzspeclib_snapshots.py"]
fn capture_fixture_snapshots() {
    let config: Value = serde_json::from_slice(
        &std::fs::read(
            std::env::var_os("MZSPECLIB_SNAPSHOT_CONFIG").expect("snapshot config required"),
        )
        .unwrap(),
    )
    .unwrap();
    let root = Path::new(config["root"].as_str().unwrap());
    let output = Path::new(config["output"].as_str().unwrap());
    for fixture in config["fixtures"].as_array().unwrap() {
        let relative = fixture.as_str().unwrap();
        eprintln!("capturing {relative}");
        let value = capture(root, relative);
        let destination = output.join(format!("{relative}.json"));
        std::fs::create_dir_all(destination.parent().unwrap()).unwrap();
        serde_json::to_writer_pretty(
            std::io::BufWriter::new(std::fs::File::create(destination).unwrap()),
            &value,
        )
        .unwrap();
    }
}
