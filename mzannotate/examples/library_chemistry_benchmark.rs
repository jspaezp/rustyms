//! Count residues and carbon atoms using serial or worker-owned reusable records.
#![allow(unused_crate_dependencies)]
use mzannotate::mzspeclib::{
    Analyte, AnalyteTarget, MzSpecLibTextParser,
    record::{MzSpecLibLibrary, MzSpecLibRecordReader, SpectrumRecord, ValueView},
};
use mzcore::{chemistry::OutputMolecularFormula, ontology::STATIC_ONTOLOGIES, prelude::*};
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

#[cfg(feature = "allocation-counting")]
#[path = "support/allocation_counting.rs"]
mod allocation_counting;
#[cfg(feature = "allocation-counting")]
#[global_allocator]
static ALLOCATOR: allocation_counting::CountingAllocator =
    allocation_counting::CountingAllocator::new();

type LegacySpectrum = mzannotate::spectrum::AnnotatedSpectrum<OutputMolecularFormula>;
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
struct Batch<'a> {
    records: Vec<SpectrumRecord<'a>>,
    active: usize,
}
impl<'a> Batch<'a> {
    fn fill<R: BufRead>(&mut self, reader: &mut MzSpecLibRecordReader<'a, R>) -> WorkResult<()> {
        self.active = 0;
        for record in &mut self.records {
            if !reader.read_into(record).map_err(|e| e.to_string())? {
                break;
            }
            self.active += 1;
        }
        Ok(())
    }
}
fn parallel<R: BufRead>(
    library: &mut MzSpecLibLibrary<'_, R>,
    workers: usize,
    batch_size: usize,
) -> WorkResult<[Counts; 2]> {
    if workers == 0 || batch_size == 0 {
        return Err("workers and batch size must be positive".into());
    }
    let mut reader = library.reader();
    thread::scope(|scope| {
        let (returned_tx, returned_rx) = mpsc::sync_channel(workers);
        let mut senders = Vec::new();
        let mut handles = Vec::new();
        for worker in 0..workers {
            let (tx, rx) = mpsc::sync_channel::<Batch<'_>>(1);
            let returned_tx = returned_tx.clone();
            handles.push(scope.spawn(move || {
                let mut counts = [Counts::default(), Counts::default()];
                while let Ok(batch) = rx.recv() {
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        batch.records[..batch.active]
                            .iter()
                            .try_for_each(|record| record_counts(record, &mut counts))
                    }))
                    .unwrap_or_else(|_| Err("worker panicked".into()));
                    if returned_tx.send((worker, batch, result)).is_err() {
                        break;
                    }
                }
                counts
            }));
            senders.push(tx);
        }
        drop(returned_tx);
        // Each worker keeps one reusable batch. Ownership returns for refill;
        // record buffers and the batch vector retain capacity. No raw-text clones.
        let mut in_flight = 0;
        let mut eof = false;
        for sender in &senders {
            let mut batch = Batch {
                records: (0..batch_size).map(|_| reader.empty_record()).collect(),
                active: 0,
            };
            batch.fill(&mut reader)?;
            eof = batch.active < batch_size;
            if batch.active == 0 {
                break;
            }
            sender.send(batch).map_err(|_| "worker disconnected")?;
            in_flight += 1;
            if eof {
                break;
            }
        }
        while in_flight != 0 {
            let (worker, mut batch, result) =
                returned_rx.recv().map_err(|_| "workers disconnected")?;
            in_flight -= 1;
            result?;
            if !eof {
                batch.fill(&mut reader)?;
                eof = batch.active < batch_size;
                if batch.active != 0 {
                    senders[worker].send(batch).map_err(|_| "worker disconnected")?;
                    in_flight += 1;
                }
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
fn legacy_counts(spectrum: &LegacySpectrum, counts: &mut [Counts; 2]) -> WorkResult<()> {
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
    counts[usize::from(decoy)].analytes(&spectrum.analytes)
}
fn fill_legacy(
    input: &mut impl Iterator<Item = WorkResult<LegacySpectrum>>,
    batch: &mut Vec<LegacySpectrum>,
    batch_size: usize,
) -> WorkResult<()> {
    for item in input.take(batch_size) {
        batch.push(item?);
    }
    Ok(())
}
fn parallel_legacy<I: Iterator<Item = WorkResult<LegacySpectrum>>>(
    mut input: I,
    workers: usize,
    batch_size: usize,
) -> WorkResult<[Counts; 2]> {
    if workers == 0 || batch_size == 0 {
        return Err("workers and batch size must be positive".into());
    }
    thread::scope(|scope| {
        let (returned_tx, returned_rx) = mpsc::sync_channel(workers);
        let mut senders = Vec::new();
        let mut handles = Vec::new();
        for worker in 0..workers {
            let (tx, rx) = mpsc::sync_channel::<Vec<LegacySpectrum>>(1);
            let returned_tx = returned_tx.clone();
            handles.push(scope.spawn(move || {
                let mut counts = [Counts::default(), Counts::default()];
                while let Ok(mut batch) = rx.recv() {
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        batch.iter().try_for_each(|record| legacy_counts(record, &mut counts))
                    }))
                    .unwrap_or_else(|_| Err("worker panicked".into()));
                    batch.clear();
                    if returned_tx.send((worker, batch, result)).is_err() {
                        break;
                    }
                }
                counts
            }));
            senders.push(tx);
        }
        drop(returned_tx);
        // Reuse the outer Vec; clearing it drops all spectrum-owned allocations.
        let mut in_flight = 0;
        let mut eof = false;
        for sender in &senders {
            let mut batch = Vec::with_capacity(batch_size);
            fill_legacy(&mut input, &mut batch, batch_size)?;
            eof = batch.len() < batch_size;
            if batch.is_empty() {
                break;
            }
            sender.send(batch).map_err(|_| "worker disconnected")?;
            in_flight += 1;
            if eof {
                break;
            }
        }
        while in_flight != 0 {
            let (worker, mut batch, result) =
                returned_rx.recv().map_err(|_| "workers disconnected")?;
            in_flight -= 1;
            result?;
            if !eof {
                fill_legacy(&mut input, &mut batch, batch_size)?;
                eof = batch.len() < batch_size;
                if !batch.is_empty() {
                    senders[worker].send(batch).map_err(|_| "worker disconnected")?;
                    in_flight += 1;
                }
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
    if !(4..=5).contains(&args.len()) {
        return Err("usage: library_chemistry_benchmark legacy|record WORKERS PATH(.gz|.txt) [BATCH=1]; WORKERS=0 means serial".into());
    }
    let workers: usize = args[2].parse()?;
    let batch_size: usize = args.get(4).map_or(Ok(1), |s| s.parse())?;
    if batch_size == 0 {
        return Err("batch size must be positive".into());
    }
    if args[1] != "record" && args[1] != "legacy" {
        return Err("unknown mode".into());
    }
    std::hint::black_box(&*STATIC_ONTOLOGIES);
    #[cfg(feature = "allocation-counting")]
    let allocation_start = ALLOCATOR.begin();
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
        let parser = MzSpecLibTextParser::open(input, None, &STATIC_ONTOLOGIES)?;
        if workers == 0 {
            for item in parser {
                legacy_counts(&item?, &mut totals)?;
            }
        } else {
            totals = parallel_legacy(
                parser.map(|item| item.map_err(|e| e.to_string())),
                workers,
                batch_size,
            )?;
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
            totals = parallel(&mut library, workers, batch_size)?;
        }
    }
    let elapsed = start.elapsed();
    #[cfg(feature = "allocation-counting")]
    let allocations = ALLOCATOR.finish(allocation_start);
    println!(
        "mode={} workers={workers} batch={batch_size} seconds={:.6}",
        args[1],
        elapsed.as_secs_f64()
    );
    #[cfg(feature = "allocation-counting")]
    println!("allocations {allocations:?}");
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
            let counts = parallel(&mut library, workers, 1).unwrap();
            assert_eq!(counts[0].carbon, 13); // A=3, C=3+2(CAM), M=5; oxidation adds no carbon.
            assert_eq!(counts[0].residues.iter().sum::<u64>(), 3);
            assert_eq!(counts[0].spectra, 1);
            assert_eq!(counts[1], Counts::default());
        }
    }

    #[test]
    fn partial_batches_do_not_recount_stale_records() {
        let (header, spectrum) = LIBRARY.split_once("<Spectrum=1>").unwrap();
        let mut input = header.to_owned();
        for key in 1..=9 {
            input.push_str(&format!("<Spectrum={key}>{spectrum}"));
        }
        for workers in [1, 2, 4] {
            for batch_size in [1, 2, 8, 32] {
                let mut library =
                    MzSpecLibLibrary::open(input.as_bytes(), None, &STATIC_ONTOLOGIES).unwrap();
                let counts = parallel(&mut library, workers, batch_size).unwrap();
                assert_eq!(counts[0].spectra, 9);
                assert_eq!(counts[0].analytes, 9);
                assert_eq!(counts[0].carbon, 9 * 13);
                assert_eq!(counts[0].residues.iter().sum::<u64>(), 9 * 3);
                assert_eq!(counts[1], Counts::default());
                let parser =
                    MzSpecLibTextParser::open(input.as_bytes(), None, &STATIC_ONTOLOGIES).unwrap();
                let legacy = parallel_legacy(
                    parser.map(|item| item.map_err(|e| e.to_string())),
                    workers,
                    batch_size,
                )
                .unwrap();
                assert_eq!(legacy, counts);
            }
        }
    }

    #[test]
    fn legacy_worker_and_producer_errors_terminate() {
        let parser =
            MzSpecLibTextParser::open(LIBRARY.as_bytes(), None, &STATIC_ONTOLOGIES).unwrap();
        let missing_origin = parser.map(|item| {
            let mut spectrum = item.map_err(|e| e.to_string())?;
            spectrum.description.params.clear();
            Ok(spectrum)
        });
        assert!(
            parallel_legacy(missing_origin, 4, 1)
                .unwrap_err()
                .contains("missing origin")
        );
        let parser =
            MzSpecLibTextParser::open(LIBRARY.as_bytes(), None, &STATIC_ONTOLOGIES).unwrap();
        let interrupted = parser
            .map(|item| item.map_err(|e| e.to_string()))
            .chain(std::iter::once(Err("input interrupted".into())));
        assert_eq!(
            parallel_legacy(interrupted, 4, 1).unwrap_err(),
            "input interrupted"
        );
    }

    #[test]
    fn worker_error_does_not_strand_other_workers() {
        let input = LIBRARY.replace("AC[UNIMOD:4]M[UNIMOD:35]/2", "???");
        let mut library =
            MzSpecLibLibrary::open(input.as_bytes(), None, &STATIC_ONTOLOGIES).unwrap();
        assert!(parallel(&mut library, 4, 8).unwrap_err().contains("ProForma"));
    }
}
