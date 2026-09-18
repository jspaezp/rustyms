//! Count residues and carbon atoms using serial or worker-owned reusable records.
#![allow(unused_crate_dependencies)]
use mzannotate::mzspeclib::{
    Analyte, AnalyteTarget, MzSpecLibTextParser,
    record::{MzSpecLibLibrary, SpectrumRecord, ValueView},
};
use mzcore::{ontology::STATIC_ONTOLOGIES, prelude::*};
use mzcv::curie;
use mzdata::params::{ControlledVocabulary, ParamValue};
use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader},
    sync::mpsc,
    thread,
    time::Instant,
};

type WorkResult<T> = Result<T, String>;
#[derive(Default, Debug, PartialEq, Eq)]
struct Counts {
    spectra: u64,
    analytes: u64,
    carbon: u64,
    residues: [u64; 26],
}
impl Counts {
    fn add(&mut self, other: &Self) {
        self.spectra += other.spectra;
        self.analytes += other.analytes;
        self.carbon += other.carbon;
        for (a, b) in self.residues.iter_mut().zip(other.residues) {
            *a += b;
        }
    }
    fn analytes(&mut self, analytes: &[Analyte]) -> WorkResult<()> {
        self.spectra += 1;
        for analyte in analytes {
            let AnalyteTarget::PeptidoformIon(ion) = &analyte.target else {
                return Err("benchmark requires peptidoform analytes".into());
            };
            self.analytes += 1;
            for peptide in ion.peptidoforms() {
                for residue in peptide.sequence() {
                    let code = residue.aminoacid.one_letter_code().ok_or("missing residue code")?;
                    if !code.is_ascii_uppercase() {
                        return Err("noncanonical residue code".into());
                    }
                    self.residues[code as usize - 'A' as usize] += 1;
                }
            }
            // Decode chemistry, including modifications and isotope-specific carbon.
            // Refuse ambiguity rather than silently selecting one candidate formula.
            let formula = analyte
                .target
                .formulas()
                .single()
                .ok_or("expected one molecular formula")?;
            for (element, _, count) in formula.elements() {
                if *element == Element::C {
                    self.carbon += u64::try_from(*count).map_err(|_| "negative carbon count")?;
                }
            }
        }
        Ok(())
    }
}
fn record_counts(record: &SpectrumRecord<'_>, counts: &mut [Counts; 2]) -> WorkResult<()> {
    if record
        .resolution_diagnostics()
        .map_err(|e| e.to_string())?
        .next()
        .is_some()
    {
        return Err("unresolved metadata".into());
    }
    let attrs = record.attributes().map_err(|e| e.to_string())?;
    let (mut predicted, mut decoy) = (false, false);
    for a in attrs.by_accession(curie!(MS:1003072)) {
        let ValueView::Term(t) = a.value().map_err(|e| e.to_string())? else {
            return Err("non-term origin".into());
        };
        match t.accession {
            Some(c) if c == curie!(MS:1003074) => predicted = true,
            Some(c) if c == curie!(MS:1003195) => decoy = true,
            _ => return Err("unknown origin".into()),
        }
    }
    if !predicted && !decoy {
        return Err("missing origin".into());
    }
    counts[usize::from(decoy)].analytes(record.analytes().map_err(|e| e.to_string())?)
}
fn parallel<R: BufRead>(
    library: &mut MzSpecLibLibrary<'_, R>,
    workers: usize,
) -> WorkResult<[Counts; 2]> {
    let mut reader = library.reader();
    thread::scope(|scope| {
        let (returned_tx, returned_rx) = mpsc::sync_channel(workers);
        let mut senders = Vec::new();
        let mut handles = Vec::new();
        for worker in 0..workers {
            let (tx, rx) = mpsc::sync_channel::<SpectrumRecord<'_>>(1);
            let returned_tx = returned_tx.clone();
            handles.push(scope.spawn(move || {
                let mut counts = [Counts::default(), Counts::default()];
                while let Ok(record) = rx.recv() {
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        record_counts(&record, &mut counts)
                    }))
                    .unwrap_or_else(|_| Err("worker panicked".into()));
                    if returned_tx.send((worker, record, result)).is_err() {
                        break;
                    }
                }
                counts
            }));
            senders.push(tx);
        }
        drop(returned_tx);
        // Exactly one record per worker. Workers return ownership for refill;
        // String/Vec capacities survive each round trip. No record clones.
        let mut in_flight = 0;
        for sender in &senders {
            let mut record = reader.empty_record();
            if !reader.read_into(&mut record).map_err(|e| e.to_string())? {
                break;
            }
            sender.send(record).map_err(|_| "worker disconnected")?;
            in_flight += 1;
        }
        let mut eof = false;
        while in_flight != 0 {
            let (worker, mut record, result) =
                returned_rx.recv().map_err(|_| "workers disconnected")?;
            in_flight -= 1;
            result?;
            if !eof && reader.read_into(&mut record).map_err(|e| e.to_string())? {
                senders[worker].send(record).map_err(|_| "worker disconnected")?;
                in_flight += 1;
            } else {
                eof = true;
            }
        }
        drop(senders);
        let mut totals = [Counts::default(), Counts::default()];
        for handle in handles {
            let counts = handle.join().map_err(|_| "worker panicked")?;
            for (total, count) in totals.iter_mut().zip(counts) {
                total.add(&count);
            }
        }
        Ok(totals)
    })
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: library_chemistry_benchmark legacy|record WORKERS PATH(.gz|.txt); WORKERS=0 means serial".into());
    }
    let workers: usize = args[2].parse()?;
    if args[1] != "record" && args[1] != "legacy" {
        return Err("unknown mode".into());
    }
    if args[1] == "legacy" && workers != 0 {
        return Err("legacy mode is serial".into());
    }
    std::hint::black_box(&*STATIC_ONTOLOGIES);
    let start = Instant::now();
    let file = BufReader::with_capacity(256 * 1024, File::open(&args[3])?);
    let input: Box<dyn BufRead> = if args[3].ends_with(".gz") {
        Box::new(BufReader::with_capacity(
            256 * 1024,
            flate2::read::MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(file)
    };
    let mut totals = [Counts::default(), Counts::default()];
    if args[1] == "legacy" {
        for item in MzSpecLibTextParser::open(input, None, &STATIC_ONTOLOGIES)? {
            let spectrum = item?;
            let (mut predicted, mut decoy) = (false, false);
            for p in &spectrum.description.params {
                if p.curie().is_some_and(|c| {
                    c.controlled_vocabulary == ControlledVocabulary::MS && c.accession == 1003072
                }) {
                    match p.value.as_str().split('|').next() {
                        Some("MS:1003074") => predicted = true,
                        Some("MS:1003195") => decoy = true,
                        _ => return Err("unknown origin".into()),
                    }
                }
            }
            if !predicted && !decoy {
                return Err("missing origin".into());
            }
            totals[usize::from(decoy)].analytes(&spectrum.analytes)?;
        }
    } else {
        let mut library = MzSpecLibLibrary::open(input, None, &STATIC_ONTOLOGIES)?;
        if workers == 0 {
            let mut reader = library.reader();
            let mut record = reader.empty_record();
            while reader.read_into(&mut record)? {
                record_counts(&record, &mut totals)?;
            }
        } else {
            totals = parallel(&mut library, workers)?;
        }
    }
    println!(
        "mode={} workers={workers} seconds={:.6}",
        args[1],
        start.elapsed().as_secs_f64()
    );
    for (label, c) in ["target", "decoy"].into_iter().zip(totals) {
        println!(
            "{label} spectra={} analytes={} residues={} carbon={} histogram={:?}",
            c.spectra,
            c.analytes,
            c.residues.iter().sum::<u64>(),
            c.carbon,
            c.residues
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const LIBRARY: &str = "<mzSpecLib>\n<AttributeSet Spectrum=all>\nMS:1003072|spectrum origin type=MS:1003074|predicted spectrum\n<Spectrum=1>\n<Analyte=1>\nMS:1003270|proforma peptidoform ion notation=AC[UNIMOD:4]M[UNIMOD:35]/2\n<Peaks>\n100\t1\t?\n";

    #[test]
    fn modified_analyte_counts_agree_across_worker_counts() {
        for workers in [1, 2, 4] {
            let mut library =
                MzSpecLibLibrary::open(LIBRARY.as_bytes(), None, &STATIC_ONTOLOGIES).unwrap();
            let counts = parallel(&mut library, workers).unwrap();
            assert_eq!(counts[0].carbon, 13); // A=3, C=3+2(CAM), M=5; oxidation adds no carbon.
            assert_eq!(counts[0].residues.iter().sum::<u64>(), 3);
            assert_eq!(counts[0].spectra, 1);
            assert_eq!(counts[1], Counts::default());
        }
    }

    #[test]
    fn worker_error_does_not_strand_other_workers() {
        let input = LIBRARY.replace("AC[UNIMOD:4]M[UNIMOD:35]/2", "???");
        let mut library =
            MzSpecLibLibrary::open(input.as_bytes(), None, &STATIC_ONTOLOGIES).unwrap();
        assert!(parallel(&mut library, 4).unwrap_err().contains("ProForma"));
    }
}
