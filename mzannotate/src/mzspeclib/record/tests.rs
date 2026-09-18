use super::*;
use mzcore::ontology::STATIC_ONTOLOGIES;
use mzcv::curie;
use mzdata::mzpeaks::prelude::PeakCollection;
use std::io::{BufReader, Cursor};

fn library(text: &str) -> MzSpecLibLibrary<'static, Cursor<Vec<u8>>> {
    MzSpecLibLibrary::open(
        Cursor::new(text.as_bytes().to_vec()),
        Some("fixture.txt".into()),
        &STATIC_ONTOLOGIES,
    )
    .unwrap()
}

#[test]
fn metadata_rejects_decoy_without_peaks_or_chemistry() {
    let mut lib = library(
        "<mzSpecLib>\n<AttributeSet Spectrum=all>\nMS:1003072|spectrum origin type=MS:1003073|observed spectrum\n<AttributeSet Spectrum=decoy>\nMS:1003072|spectrum origin type=MS:1003195|shuffle-and-reposition decoy spectrum\nMS:1003072|spectrum origin type=MS:1003074|predicted spectrum\n<Spectrum=7>\nMS:1003212|library attribute set name=decoy\nMS:1003212|library attribute set name=undefined\nMS:1000894|retention time=bad\n<Analyte=1>\nMS:1003270|proforma peptidoform ion notation=???\n<Peaks>\ninvalid\n",
    );
    let mut reader = lib.reader();
    let mut record = reader.empty_record();
    assert!(reader.read_into(&mut record).unwrap());
    let attrs = record.attributes().unwrap();
    let origins: Vec<_> = attrs.by_accession(curie!(MS:1003072)).collect();
    assert_eq!(origins.len(), 2);
    assert!(
        matches!(origins[0].value().unwrap(), ValueView::Term(t) if t.accession == Some(curie!(MS:1003195)))
    );
    assert_eq!(origins[0].origin().attribute_set, Some("decoy"));
    assert!(attrs.by_accession(curie!(MS:1000894)).next().unwrap().to_f64().is_err());
    assert_eq!(record.resolution_diagnostics().unwrap().count(), 1);
    assert!(record.analytes.completed.get().is_none());
    assert!(record.peaks.completed.get().is_none());
    assert!(record.annotations.completed.get().is_none());
    assert!(record.analytes().is_err());
    assert!(record.peaks().is_err());
    assert_eq!(
        record.analyte_scopes().next().unwrap().local_attributes().len(),
        1
    );
    assert_eq!(record.key(), Some(7));
}

#[test]
fn peaks_preserve_fields_and_temporary_view_lifetimes() {
    let mut lib = library(
        "<mzSpecLib>\r\n<Spectrum=1>\r\n<Peaks>\r\n100\t1\r\n101\t2\t\r\n102\t3\t?\t\r\n103\t4\t y2 \t1,2\textra\t\r\n",
    );
    let record = lib.reader().records().next().unwrap().unwrap();
    let first = record.peaks().unwrap().get(0).unwrap();
    assert_eq!(first.annotation_field(), SuppliedField::Absent);
    assert_eq!(
        record.peaks().unwrap().get(1).unwrap().annotation_field(),
        SuppliedField::Present("")
    );
    assert_eq!(
        record.peaks().unwrap().get(2).unwrap().annotation_field(),
        SuppliedField::Present("?")
    );
    let last = record.peaks().unwrap().get(3).unwrap();
    assert_eq!(last.annotation_field(), SuppliedField::Present(" y2 "));
    assert_eq!(last.aggregation_field(), SuppliedField::Present("1,2"));
    assert_eq!(last.extra_columns().collect::<Vec<_>>(), ["extra", ""]);
    assert_eq!(last.source_row(), 3);
    assert_eq!(last.source_position().line, 6);
    let peaks = record.peaks().unwrap();
    assert_eq!(peaks.into_iter().count(), 4);
    assert_eq!(peaks.mz().as_ptr(), record.peaks().unwrap().mz().as_ptr());
}

#[test]
fn reports_keep_syntax_when_chemistry_fails_and_all_field_errors() {
    let mut lib = library(
        "<mzSpecLib>\n<Spectrum=1>\n<Analyte=1>\nMS:1003270|proforma peptidoform ion notation=???\n<Peaks>\n100\t1\ty2^2\n101\t2\tnot_an_annotation\n102\t3\t?\n",
    );
    let record = lib.reader().records().next().unwrap().unwrap();
    let row = record.resolved_annotations().unwrap().get(0).unwrap();
    assert_eq!(row.iter().len(), 1);
    assert!(row.iter().next().unwrap().fragment.is_none());
    assert_eq!(row.annotation_field(), SuppliedField::Present("y2^2"));
    let report = record.resolved_annotations().unwrap();
    assert_eq!(report.len(), 3);
    assert_eq!(report.diagnostics().count(), 2);
    assert!(report.diagnostics().any(|d| d.source_row == 1 && d.alternative.is_none()));
    assert!(report.diagnostics().any(|d| d.error.kind() == RecordErrorKind::Unresolved));
}

#[test]
fn cache_reuses_capacity_on_success_error_skipped_success() {
    let mut slot = Reusable::<Vec<u64>>::default();
    let pointer = slot
        .get(|v| {
            v.extend(0..64);
            Ok(())
        })
        .unwrap()
        .as_ptr();
    slot.reset();
    let source = SourcePosition {
        source: SourceId::new(),
        line: 0,
        byte_offset: 0,
    };
    assert!(
        slot.get(|v| {
            v.extend(0..32);
            Err(RecordError::new(
                RecordErrorKind::Malformed,
                "expected",
                source,
            ))
        })
        .is_err()
    );
    assert_eq!(slot.completed.get().unwrap().storage.as_ptr(), pointer);
    assert!(slot.get(|_| panic!("cached failure must not retry")).is_err());
    slot.reset();
    slot.reset();
    assert_eq!(
        slot.get(|v| {
            v.extend(0..16);
            Ok(())
        })
        .unwrap()
        .as_ptr(),
        pointer
    );
}

#[test]
fn reuse_changes_generation_not_owner_and_invalidates_groups() {
    let text = "<mzSpecLib>\n<Spectrum=1>\n[1]MS:1000894|retention time=0\n[1]UO:0000000|unit=UO:0000010|second\n<Peaks>\n100\t1\n<Spectrum=2>\n[1]MS:1000894|retention time=2\n<Peaks>\n100\tbad\n<Spectrum=3>\n<Peaks>\n100\t3\n";
    let mut lib = library(text);
    let mut reader = lib.reader();
    let mut record = reader.empty_record();
    reader.read_into(&mut record).unwrap();
    let owner = record.owner;
    let generation = record.generation;
    let id = record.attributes().unwrap().iter().next().unwrap().group_id().unwrap();
    let pointer = record.peaks().unwrap().mz().as_ptr();
    let raw_pointer = record.raw.text.as_ptr();
    reader.read_into(&mut record).unwrap();
    assert_eq!(record.owner, owner);
    assert_ne!(record.generation, generation);
    assert_eq!(record.raw.text.as_ptr(), raw_pointer);
    assert!(record.attributes().unwrap().group(id).is_none());
    assert!(record.peaks().is_err());
    reader.read_into(&mut record).unwrap();
    assert_eq!(record.peaks().unwrap().mz().as_ptr(), pointer);
    assert!(!reader.read_into(&mut record).unwrap());
    assert!(!record.is_loaded());
    assert_eq!(record.key(), None);
}

#[test]
fn retained_records_and_foreign_group_ids() {
    let mut lib = library(
        "<mzSpecLib>\n<Spectrum=1>\n[2]MS:1000894|retention time=0\n<Peaks>\n<Spectrum=2>\n[2]MS:1000894|retention time=1\n<Peaks>",
    );
    let (first, second) = {
        let mut reader = lib.reader();
        let mut records = reader.records();
        (
            records.next().unwrap().unwrap(),
            records.next().unwrap().unwrap(),
        )
    };
    let id = first.attributes().unwrap().iter().next().unwrap().group_id().unwrap();
    assert!(second.attributes().unwrap().group(id).is_none());
    assert_eq!(first.key(), Some(1));
    assert_eq!(second.key(), Some(2));
}

#[test]
fn raw_value_errors_do_not_poison_other_occurrences() {
    let mut lib = library(
        "<mzSpecLib>\n<Spectrum=1>\nMS:1000894|retention time=\"unclosed\nMS:1003072|spectrum origin type=MS:9999999|unknown origin\n<Peaks>\n",
    );
    let record = lib.reader().records().next().unwrap().unwrap();
    let attrs = record.attributes().unwrap();
    assert!(attrs.by_accession(curie!(MS:1000894)).next().unwrap().value().is_err());
    assert!(
        matches!(attrs.by_accession(curie!(MS:1003072)).next().unwrap().value().unwrap(), ValueView::Term(t) if t.accession == Some(curie!(MS:9999999)))
    );
}

#[test]
fn tiny_buffers_recovery_and_grouped_claims() {
    let text = "<mzSpecLib>\n<AttributeSet Spectrum=unit>\nUO:0000000|unit=UO:0000010|second\n<Spectrum=bad>\n<Peaks>\n<Spectrum=2>\n[8]MS:1003212|library attribute set name=unit\n[8]MS:1000894|retention time=0\n<Peaks>\n";
    let mut lib = MzSpecLibLibrary::open(
        BufReader::with_capacity(1, text.as_bytes()),
        None,
        &STATIC_ONTOLOGIES,
    )
    .unwrap();
    let mut reader = lib.reader();
    let mut record = reader.empty_record();
    assert!(reader.read_into(&mut record).is_err());
    assert!(reader.read_into(&mut record).unwrap());
    let attrs = record.attributes().unwrap();
    let rt = attrs.by_accession(curie!(MS:1000894)).next().unwrap();
    assert_eq!(rt.to_f64().unwrap(), 0.0);
    assert_eq!(
        rt.group().unwrap().by_accession(curie!(UO:0000000)).count(),
        1
    );
    assert!(!reader.read_into(&mut record).unwrap());
}

#[test]
fn all_shipped_data_records_load_and_numeric_peaks_decode() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../data");
    for entry in std::fs::read_dir(root).unwrap() {
        let path = entry.unwrap().path();
        if !path.to_string_lossy().to_ascii_lowercase().ends_with(".mzspeclib.txt") {
            continue;
        }
        let input = BufReader::new(std::fs::File::open(&path).unwrap());
        let mut lib =
            MzSpecLibLibrary::open(input, Some(path.clone()), &STATIC_ONTOLOGIES).unwrap();
        for record in lib.reader().records() {
            let record = record.unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            record.attributes().unwrap();
            record.peaks().unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        }
    }
}

#[test]
fn modern_materialization_uses_decoded_sections() {
    let mut lib = library(
        "<mzSpecLib>\n<Spectrum=1>\nMS:1003061|library spectrum name=test\n<Analyte=1>\nMS:1003270|proforma peptidoform ion notation=PEPTIDE/2\n[1]MS:1000885|protein accession=P1\nMS:1000888|stripped peptide sequence=PEPTIDE\n<Interpretation=1>\nMS:1003163|analyte mixture members=1\nMS:1002357|PSM-level probability=0.95\n<Peaks>\n100\t1\ty2\n",
    );
    let record = lib.reader().records().next().unwrap().unwrap();
    let spectrum = record.materialize().unwrap();
    assert_eq!(spectrum.key, 1);
    assert_eq!(spectrum.description.id, "test");
    assert_eq!(spectrum.analytes.len(), 1);
    assert_eq!(spectrum.analytes[0].proteins[0].accession.as_deref(), Some("P1"));
    assert!(spectrum.analytes[0].params.iter().any(|p| p.name == "stripped peptide sequence"));
    assert_eq!(spectrum.interpretations[0].probability, Some(0.95));
    assert_eq!(spectrum.interpretations[0].analyte_refs, [1]);
    assert_eq!(spectrum.peaks.len(), 1);
    let consumed = record.into_annotated().unwrap();
    assert_eq!(consumed.peaks, spectrum.peaks);
    assert_eq!(consumed.analytes, spectrum.analytes);
    assert_eq!(consumed.interpretations, spectrum.interpretations);
}

#[test]
fn indexed_reads_restore_position_and_work_after_eof() {
    let mut lib = library(
        "<mzSpecLib>\n<Spectrum=3>\nMS:1003061|library spectrum name=first\nMS:1003057|scan number=42\n<Peaks>\n100\t1\n<Spectrum=9>\nMS:1003061|library spectrum name=second\n<Peaks>\n101\t2\n",
    );
    let mut reader = lib.reader();
    let mut record = reader.empty_record();
    assert_eq!(reader.build_index().unwrap(), 2);
    assert!(reader.read_into(&mut record).unwrap());
    assert_eq!(record.key(), Some(3));
    assert!(reader.read_into(&mut record).unwrap());
    assert_eq!(record.key(), Some(9));
    assert!(!reader.read_into(&mut record).unwrap());
    assert!(reader.read_by_scan_number_into(42, &mut record).unwrap());
    assert_eq!(record.key(), Some(3));
    assert_eq!(record.source_position().unwrap().line, 1);
    assert!(reader.read_by_name_into("second", &mut record).unwrap());
    assert_eq!(record.key(), Some(9));
    assert!(reader.read_by_key_into(3, &mut record).unwrap());
    assert_eq!(record.peaks().unwrap().mz(), [100.0]);
    assert!(!reader.read_by_index_into(2, &mut record).unwrap());
    assert!(!record.is_loaded());
}

#[test]
fn source_ids_are_opaque_and_distinct_between_libraries() {
    let first = library("<mzSpecLib>\n<Spectrum=1>\n<Peaks>\n");
    let second = library("<mzSpecLib>\n<Spectrum=1>\n<Peaks>\n");
    assert_ne!(first.metadata.source, second.metadata.source);
}

#[test]
fn compound_losses_and_isotopes_are_not_singletons() {
    let losses =
        crate::fragment::PeakAnnotation::parse_all("y2-H2O-NH3", &STATIC_ONTOLOGIES).unwrap();
    assert_eq!(losses[0].neutral_losses().len(), 2);
    let isotopes =
        crate::fragment::PeakAnnotation::parse_all("y12+2i13C+i15N", &STATIC_ONTOLOGIES).unwrap();
    assert_eq!(isotopes[0].isotopes().len(), 2);
    assert_eq!(isotopes[0].isotopes()[0].0, 2);
}

#[test]
fn interpreted_properties_use_inheritance_and_retain_reported_metadata() {
    let mut lib = library(
        "<mzSpecLib>\n<AttributeSet Analyte=all>\nMS:1003270|proforma peptidoform ion notation=PEPTIDE/2\n<AttributeSet Interpretation=all>\nMS:1002357|PSM-level probability=0.9\n<Spectrum=1>\n<Analyte=1>\nMS:1003053|theoretical monoisotopic m/z=123.4\n<Interpretation=1>\n<Peaks>\n",
    );
    let record = lib.reader().records().next().unwrap().unwrap();
    assert!(matches!(
        record.analytes().unwrap()[0].target,
        crate::mzspeclib::AnalyteTarget::PeptidoformIon(_)
    ));
    assert_eq!(record.interpretations().unwrap()[0].probability, Some(0.9));
    assert_eq!(
        record
            .analyte_scopes()
            .next()
            .unwrap()
            .attributes()
            .unwrap()
            .by_accession(curie!(MS:1003053))
            .next()
            .unwrap()
            .raw_value(),
        "123.4"
    );
}

#[test]
fn records_move_to_scoped_workers_and_return_for_reuse() {
    fn send<T: Send>() {}
    fn sync<T: Sync>() {}
    send::<SpectrumRecord<'_>>();
    sync::<LibraryContext<'_>>();
    sync::<LibraryMetadata>();

    let mut lib = library(
        "<mzSpecLib>\n<AttributeSet Spectrum=all>\nMS:1003072|spectrum origin type=MS:1003074|predicted spectrum\nMS:1000894|retention time=\"bad\n<Spectrum=1>\n<Peaks>\n100\t1\n<Spectrum=2>\n<Peaks>\n200\t2\n<Spectrum=3>\n<Peaks>\n300\t3\n",
    );
    // Header values are decoded at open; malformed values remain per-occurrence errors.
    let mut reader = lib.reader();
    let mut first = reader.empty_record();
    let mut second = reader.empty_record();
    assert!(reader.read_into(&mut first).unwrap());
    assert!(reader.read_into(&mut second).unwrap());
    let mut returned = std::thread::scope(|scope| {
        let handles: Vec<_> = [first, second]
            .into_iter()
            .map(|record| {
                scope.spawn(move || {
                    let attrs = record.attributes().unwrap();
                    assert!(
                        attrs.by_accession(curie!(MS:1000894)).next().unwrap().value().is_err()
                    );
                    assert!(matches!(
                        attrs.by_accession(curie!(MS:1003072)).next().unwrap().value(),
                        Ok(ValueView::Term(_))
                    ));
                    record.peaks().unwrap();
                    record
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect::<Vec<_>>()
    });
    let record = &mut returned[0];
    let raw_pointer = record.raw_text().as_ptr();
    let peak_pointer = record.peaks().unwrap().mz().as_ptr();
    assert!(reader.read_into(record).unwrap());
    assert_eq!(record.key(), Some(3));
    assert_eq!(record.raw_text().as_ptr(), raw_pointer);
    assert_eq!(record.peaks().unwrap().mz().as_ptr(), peak_pointer);
    assert_eq!(record.peaks().unwrap().mz(), [300.0]);
}
