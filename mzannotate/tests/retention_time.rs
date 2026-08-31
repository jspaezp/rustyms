//! `MS:1000894|retention time` and `MS:1000896|normalized retention time` are
//! normalised into `mzdata`'s `ScanEvent::start_time`, which is documented as
//! minutes (`mzdata/src/spectrum/scan_properties.rs`). Every path has to land
//! on that unit, and a retention time has to survive whatever else its writer
//! grouped alongside it.
use mzannotate::mzspeclib::MzSpecLibTextParser;

/// `1719.50946` s and `28.658491` min are the same instant. Spectrum 3 carries
/// the same instant again with two unrelated siblings in the group, which is
/// how SpectraST writes a consensus retention time.
const LIBRARY: &str = "<mzSpecLib>
MS:1003186|library format version=1.0
<Spectrum=1>
[1]MS:1000894|retention time=1719.509460
[1]UO:0000000|unit=UO:0000010|second
<Peaks>

<Spectrum=2>
[1]MS:1000894|retention time=28.658491
[1]UO:0000000|unit=UO:0000031|minute
<Peaks>

<Spectrum=3>
[1]MS:1000894|retention time=1719.509460
[1]UO:0000000|unit=UO:0000010|second
[1]MS:1003174|attribute maximum=1800.0
[1]MS:1003175|attribute minimum=1600.0
<Peaks>
";

/// `28.658491` minutes, the value all three spectra state.
const EXPECTED_MINUTES: f64 = 28.658_491;

fn start_times() -> Vec<f64> {
    MzSpecLibTextParser::open(LIBRARY.as_bytes(), None, &mzcore::ontology::STATIC_ONTOLOGIES)
        .unwrap()
        .map(|s| s.unwrap().description.acquisition.scans[0].start_time)
        .collect()
}

#[test]
fn every_declared_unit_lands_on_minutes() {
    let times = start_times();
    assert_eq!(times.len(), 3);
    for (i, time) in times.iter().enumerate() {
        assert!(
            (time - EXPECTED_MINUTES).abs() < 1e-4,
            "spectrum {} should be {EXPECTED_MINUTES} min, got {time}",
            i + 1
        );
    }
}

/// A retention time grouped with unrelated attributes was previously dropped
/// entirely, so whether a spectrum kept its retention time depended on what
/// else its writer put in the group.
#[test]
fn a_retention_time_survives_unrelated_group_members() {
    assert_ne!(start_times()[2], 0.0);
}
