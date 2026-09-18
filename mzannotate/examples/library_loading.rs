//! Consumer-owned decoy filtering before optional section decoding.
#![allow(unused_crate_dependencies)]
use std::{error::Error, fs::File, io::BufReader, path::PathBuf};

use mzannotate::mzspeclib::record::{MzSpecLibLibrary, ValueView};
use mzcore::ontology::STATIC_ONTOLOGIES;
use mzcv::curie;

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(
        args.next()
            .ok_or("usage: library_loading PATH [--peaks] [--annotations]")?,
    );
    let flags: Vec<_> = args.collect();
    if flags.iter().any(|s| s != "--peaks" && s != "--annotations") {
        return Err("unknown option".into());
    }
    let annotations = flags.iter().any(|s| s == "--annotations");
    let peaks = annotations || flags.iter().any(|s| s == "--peaks");
    let mut library = MzSpecLibLibrary::open(
        BufReader::new(File::open(&path)?),
        Some(path),
        &STATIC_ONTOLOGIES,
    )?;
    let mut reader = library.reader();
    let mut record = reader.empty_record();
    let (mut retained, mut rejected, mut decoded_peaks, mut alternatives) = (0, 0, 0, 0);
    while reader.read_into(&mut record)? {
        let attrs = record.attributes()?;
        let mut known_decoy = false;
        for occurrence in attrs.by_accession(curie!(MS:1003072)) {
            // Deliberately limited consumer policy for the shipped SpectraST fixture.
            // Production consumers should use their ontology ancestry classifier.
            if let Ok(ValueView::Term(term)) = occurrence.value() {
                known_decoy |= term.accession == Some(curie!(MS:1003195));
            }
        }
        for diagnostic in record.resolution_diagnostics()? {
            eprintln!("key={:?}: {diagnostic}", record.key());
        }
        if known_decoy {
            rejected += 1;
            continue; // Neither peaks() nor analytes() was requested.
        }
        retained += 1; // Unclassified, NOT an assertion that this is a target.
        let names: Vec<_> = attrs.by_accession(curie!(MS:1003061)).map(|v| v.raw_value()).collect();
        println!(
            "key={:?} names={names:?} line={}",
            record.key(),
            record.source_position().unwrap().line + 1
        );
        if peaks {
            decoded_peaks += record.peaks()?.iter().count();
        }
        if annotations {
            let report = record.resolved_annotations()?;
            for diagnostic in report.diagnostics() {
                eprintln!(
                    "key={:?} peak={}: {}",
                    record.key(),
                    diagnostic.source_row,
                    diagnostic.error
                );
            }
            alternatives += report.iter().map(|row| row.iter().len()).sum::<usize>();
        }
    }
    println!(
        "retained={retained} rejected_known_decoys={rejected} decoded_peaks={decoded_peaks} reported_alternatives={alternatives}"
    );
    Ok(())
}
