# mzSpecLib aggregate benchmark

Workload: sum all peak intensities, mean reported fragment m/z over all peak rows,
and mean selected-ion precursor m/z over all spectra, separately for targets and
decoys. Means are unweighted. Every spectrum and every peak is included.

The executable is [library_benchmark.rs](../mzannotate/examples/library_benchmark.rs).
It compares the existing `MzSpecLibTextParser` with `MzSpecLibLibrary` using one
reusable record. The record API requests effective metadata and numeric peaks;
it does not request analyte chemistry or annotation resolution. The legacy API
materializes those sections as part of reading. This measures the requested
consumer workload, not equal-work full chemical decoding throughput.

This file declares predicted origins (`MS:1003074`) for targets and explicit decoy
origins (`MS:1003195`) through a named attribute set. Classification compares
accessions, not supplied CV labels or protein names. The sample policy rejects
missing/unknown origins, and gives decoy evidence precedence if both declarations
are present. This is a policy for this fixture, not a general ontology classifier.
Precursor coordinates use only `MS:1000744` (selected ion m/z). The record path
requires exactly one; neither theoretical m/z nor fragment-derived mass substitutes
for it. Intensities use each API's `f32` values accumulated as `f64`; both m/z sums
use `f64`. Counts and printed aggregates are checked across every run.

## Reproduce

From the repository root (the input is external and is not committed):

```sh
mkdir -p .local-design
gzip -dc ~/fasta/hela_gt20peps.mzspeclib.txt.gz > .local-design/hela_gt20peps.mzspeclib.txt
cargo build --release --offline --locked -p mzannotate --example library_benchmark
for mode in record legacy legacy record record legacy; do
  target/release/examples/library_benchmark "$mode" .local-design/hela_gt20peps.mzspeclib.txt
done
rm .local-design/hela_gt20peps.mzspeclib.txt
```

Build time, process startup, ontology initialization, decompression and result
printing are outside the reported timers. Opening/reading the decompressed file,
header parsing, record parsing, allocation, classification, aggregation and reader
cleanup are inside. Both use a 256 KiB input buffer. Runs are sequential, with warm
filesystem cache; this is not a cold-storage benchmark. No allocator instrumentation
or peak-memory measurement was performed.

## Measured results (2026-09-17)

macOS 26.5.1, arm64; rustc 1.95.0 (59807616e 2026-04-14). Parser
implementation: `54012002`; release profile: opt-level 3, LTO, one codegen unit.
Input: `hela_gt20peps.mzspeclib.txt.gz`, 266,335,688 compressed bytes,
1,359,658,233 decompressed bytes. Compressed SHA-256:
`d58e749ad193fcd97dda63734b18b0fc6edd4e77499f0dc1b734721449d8d636`.

Three passes per API, in order: record, legacy, legacy, record, record, legacy.

| API | Three parse + aggregate times (s) | Median (s) |
| --- | --- | --- |
| Legacy | 22.982024, 22.980978, 22.962760 | 22.980978 |
| Reusable record | 6.320962, 6.206107, 6.217299 | 6.217299 |

**3.70× speedup** for this workload;
72.9% less measured time. Both APIs completed without
errors: 948,957 spectra and 27,287,703 peaks. All six runs produced identical
counts and aggregates to the 12 decimal places printed below.

```text
target spectra=474630 peaks=13644852 intensity_sum=2415009.778824022040 fragment_mz_mean=812.012347993794 precursor_mz_mean=693.016478019259
decoy spectra=474327 peaks=13642851 intensity_sum=2415411.354313865304 fragment_mz_mean=812.422592525344 precursor_mz_mean=693.255321040978
```

A separate one-pass gzip decompression into the temporary text file took 4.164 s
(Python gzip + 1 MiB copying buffers, including output writes). This was preparation,
concurrent with the release build, not a controlled decompression benchmark. It is
excluded from every parse timing above. Streaming gzip end-to-end time was not
measured and cannot be inferred by simply adding these timings.
