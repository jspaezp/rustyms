//! Compare legacy materialization and reusable numeric record access on a text library.
#![allow(unused_crate_dependencies)]
use mzannotate::mzspeclib::{
    MzSpecLibTextParser,
    record::{MzSpecLibLibrary, ValueView},
};
use mzcore::ontology::STATIC_ONTOLOGIES;
use mzcv::curie;
use mzdata::params::{ControlledVocabulary, ParamValue};
use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader},
    time::Instant,
};

#[derive(Default, Debug)]
struct Totals {
    spectra: u64,
    peaks: u64,
    intensity: f64,
    fragment_mz: f64,
    precursor_mz: f64,
}
impl Totals {
    fn precursor(&mut self, mz: f64) {
        assert!(mz.is_finite() && mz > 0.0);
        self.spectra += 1;
        self.precursor_mz += mz;
    }
    fn peak(&mut self, mz: f64, intensity: f32) {
        assert!(mz.is_finite() && intensity.is_finite());
        self.peaks += 1;
        self.fragment_mz += mz;
        self.intensity += f64::from(intensity);
    }
}
// Fixture policy: predicted-only is target; explicit MS:1003195 is decoy.
// Fail on missing/unknown origins instead of silently counting them as targets.
#[allow(single_use_lifetimes)]
fn classify<'a>(
    origins: impl Iterator<Item = std::borrow::Cow<'a, str>>,
) -> Result<usize, Box<dyn Error>> {
    let (mut predicted, mut decoy) = (false, false);
    for value in origins {
        match value.split('|').next().unwrap_or("") {
            "MS:1003074" => predicted = true,
            "MS:1003195" => decoy = true,
            other => return Err(format!("unsupported origin {other}").into()),
        }
    }
    if decoy {
        Ok(1)
    } else if predicted {
        Ok(0)
    } else {
        Err("missing origin".into())
    }
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    let [_, mode, path] = args.as_slice() else {
        return Err("usage: library_benchmark legacy|record PATH (.txt or .gz)".into());
    };
    // Force ontology initialization outside the measured region for both paths.
    std::hint::black_box(&*STATIC_ONTOLOGIES);
    let start = Instant::now();
    let compressed = path.ends_with(".gz");
    let file = BufReader::with_capacity(256 * 1024, File::open(path)?);
    let input: Box<dyn BufRead> = if compressed {
        Box::new(BufReader::with_capacity(
            256 * 1024,
            flate2::read::MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(file)
    };
    let mut totals = [Totals::default(), Totals::default()];
    match mode.as_str() {
        "legacy" => {
            for item in MzSpecLibTextParser::open(input, None, &STATIC_ONTOLOGIES)? {
                let spectrum = item?;
                let bucket = classify(
                    spectrum
                        .description
                        .params
                        .iter()
                        .filter(|p| {
                            p.curie().is_some_and(|c| {
                                c.controlled_vocabulary == ControlledVocabulary::MS
                                    && c.accession == 1003072
                            })
                        })
                        .map(|p| p.value.as_str()),
                )?;
                let precursor =
                    spectrum.description.precursor.first().ok_or("missing precursor")?;
                if precursor.ions.len() != 1 {
                    return Err("expected one selected ion".into());
                }
                totals[bucket].precursor(precursor.ions[0].mz);
                for peak in spectrum.peaks.iter() {
                    totals[bucket].peak(peak.mz.value, peak.intensity);
                }
            }
        }
        "record" => {
            let mut library = MzSpecLibLibrary::open(input, None, &STATIC_ONTOLOGIES)?;
            let mut reader = library.reader();
            let mut record = reader.empty_record();
            while reader.read_into(&mut record)? {
                if record.resolution_diagnostics()?.next().is_some() {
                    return Err("unresolved metadata".into());
                }
                let attrs = record.attributes()?;
                let (mut predicted, mut decoy) = (false, false);
                for origin in attrs.by_accession(curie!(MS:1003072)) {
                    let ValueView::Term(term) = origin.value().map_err(Clone::clone)? else {
                        return Err("non-term origin".into());
                    };
                    match term.accession {
                        Some(c) if c == curie!(MS:1003074) => predicted = true,
                        Some(c) if c == curie!(MS:1003195) => decoy = true,
                        _ => return Err("unsupported origin".into()),
                    }
                }
                let bucket = if decoy {
                    1
                } else if predicted {
                    0
                } else {
                    return Err("missing origin".into());
                };
                let mut precursors = attrs.by_accession(curie!(MS:1000744));
                let precursor = precursors.next().ok_or("missing selected-ion m/z")?.to_f64()?;
                if precursors.next().is_some() {
                    return Err("multiple selected-ion m/z values".into());
                }
                totals[bucket].precursor(precursor);
                for peak in record.peaks()?.iter() {
                    totals[bucket].peak(peak.mz(), peak.intensity());
                }
            }
        }
        _ => return Err("mode must be legacy or record".into()),
    }
    let elapsed = start.elapsed().as_secs_f64();
    println!("mode={mode} gzip={compressed} seconds={elapsed:.6}");
    for (label, t) in ["target", "decoy"].into_iter().zip(totals) {
        println!(
            "{label} spectra={} peaks={} intensity_sum={:.12} fragment_mz_mean={:.12} precursor_mz_mean={:.12}",
            t.spectra,
            t.peaks,
            t.intensity,
            t.fragment_mz / t.peaks as f64,
            t.precursor_mz / t.spectra as f64
        );
    }
    Ok(())
}
