# Chemistry counting and scoped workers

The [reuse correction](#reuse-correction) below supersedes the original allocation
and timing results. The older measurements are retained as explicit before/after
evidence.

The record API now decodes header values once at library opening and stores them
immutably, including per-occurrence errors. Its borrowed library context is `Sync`;
`SpectrumRecord` is `Send` but remains non-`Sync`. A compile-time trait check and a
scoped-thread test verify moving records to workers and back, preserving raw-text
and numeric-buffer allocations on refill. No unsafe code or shared mutable record
cache is involved.

[library_chemistry_benchmark.rs](../mzannotate/examples/library_chemistry_benchmark.rs)
compares serial record access, worker-owned records and the branch legacy reader.
Every spectrum's **all analytes** are decoded. It counts each amino-acid letter and
carbon atoms in the resulting molecular formula, including modifications and
isotope-labelled carbon. Charge-state records count separately; this is not a
unique-peptide count. The example requires peptidoform analytes and a single formula
per analyte, failing rather than choosing among ambiguous formulas. Targets/decoys
use the same explicit origin policy as the numeric benchmark.

Each worker owns a reusable batch (default: one record). A single producer
decompresses gzip and frames records, then moves each loaded batch to a worker.
The worker resolves metadata, decodes analytes, counts chemistry and returns the
batch for refill. Counters stay on each worker and merge only at `join`; the return
channel recycles buffers, not partial sums. Storage is bounded at workers × batch
size records; retained byte capacity depends on record sizes. Larger batches reduce
channel handoffs while retaining more reusable buffers.
Channels are bounded; the library outlives scoped workers. The input reader stays
on the producer thread and does not need to implement Send. No raw buffers or decoded
analytes are cloned for dispatch. Record caches still use `OnceCell`/`RefCell`; only
the shared header is immutable. The common chemistry path now reuses nested storage and formula outputs; general
chemistry paths may still allocate (see the reuse correction below).
A processing error terminates the run without reporting partial totals as success.

Both record paths skip numeric peak decoding and fragment annotation resolution.
Legacy materialization performs those operations even though this workload does
not use them. Parallel-vs-serial record timings therefore measure the same requested
work; legacy comparisons additionally reflect avoided materialization work.

## Run

```sh
cargo build --release --offline --locked -p mzannotate --example library_chemistry_benchmark
# Direct serial processing (no worker dispatch).
target/release/examples/library_chemistry_benchmark record 0 ~/fasta/hela_gt20peps.mzspeclib.txt.gz
# Two/four worker threads, plus the reader thread.
target/release/examples/library_chemistry_benchmark record 2 ~/fasta/hela_gt20peps.mzspeclib.txt.gz
target/release/examples/library_chemistry_benchmark record 4 ~/fasta/hela_gt20peps.mzspeclib.txt.gz
# Batch 32 records per handoff, keeping at most 128 reusable records.
target/release/examples/library_chemistry_benchmark record 4 ~/fasta/hela_gt20peps.mzspeclib.txt.gz 32
# Existing materializing reader on this feature branch.
target/release/examples/library_chemistry_benchmark legacy 0 ~/fasta/hela_gt20peps.mzspeclib.txt.gz
```

Timers include gzip decompression, input reads, header/record parsing, chemical
formula calculation, thread/channel setup, result merging and reader cleanup.
Build/process startup, shared ontology initialization and printing are excluded.
The gzip/input buffers are 256 KiB. Runs are sequential with warm file cache.
`histogram` in stdout has 26 positions, A through Z. Sums are integer and therefore
independent of worker completion order.

Validation includes hand-counted modified `AC[UNIMOD:4]M[UNIMOD:35]/2` (3 residues,
13 carbons) and `PEPTIDE/2` (7 residues, 34 carbons), matching serial/1/2/4-worker
results. Example tests cover the modified count, partial final batches after buffer reuse,
and a failing worker without stranding other workers. The package's record test checks header errors remain
independent and records retain allocations across a worker round trip.

## Actual main baseline

`origin/main` was checked again on 2026-09-18 and remains
`8524b24ebfe62fa308532e61a8a12c4b3ae522d5`. Its isolated executable uses the same
`Counts::analytes` and legacy loop, with unavailable record functions removed;
parser code is unchanged. As in the numeric benchmark, main loses inherited decoy
origins. Its combined chemical counts agree, but its target/decoy split is wrong.
Main is shown separately, not used as an equivalent-correctness speedup baseline.

To build that isolated baseline (run from this repository root):

```sh
chemistry_main_dir=$(mktemp -d)
git worktree add --detach "$chemistry_main_dir" 8524b24ebfe62fa308532e61a8a12c4b3ae522d5
cp Cargo.lock "$chemistry_main_dir/Cargo.lock"
uv run --no-project python - "$chemistry_main_dir" <<'PY'
from pathlib import Path
import sys, subprocess
root = Path(sys.argv[1])
# Pin the serial-baseline harness used for the timings below.
s = subprocess.check_output(['git', 'show', 'e1aa005c:mzannotate/examples/library_chemistry_benchmark.rs'], text=True).split('#[cfg(test)]')[0]
s = s.replace('    record::{MzSpecLibLibrary, MzSpecLibRecordReader, SpectrumRecord, ValueView},\n', '').replace('use mzcv::curie;\n', '')
s = s.replace('#![allow(unused_crate_dependencies)]', '#![allow(unused_crate_dependencies, unused_imports, dead_code)]')
a, b = s.index('fn record_counts('), s.index('fn main()')
s = s[:a] + s[b:]
a = s.index('    } else {\n        let mut library = MzSpecLibLibrary::open')
b = s.index('    println!(\n        "mode={}', a)
s = s[:a] + '    }\n' + s[b:]
p = root / 'mzannotate/examples/library_chemistry_benchmark_main.rs'
p.parent.mkdir(exist_ok=True)
p.write_text(s)
p = root / 'mzannotate/Cargo.toml'
p.write_text(p.read_text() + '\n[dev-dependencies]\nflate2 = { workspace = true }\n')
PY
CARGO_TARGET_DIR="$PWD/target" cargo build \
  --manifest-path "$chemistry_main_dir/Cargo.toml" --release --offline \
  -p mzannotate --example library_chemistry_benchmark_main
target/release/examples/library_chemistry_benchmark_main legacy 0 ~/fasta/hela_gt20peps.mzspeclib.txt.gz
git worktree remove --force "$chemistry_main_dir"
```

Main has no committed lockfile; its build starts with this workspace's lockfile and
resolves it offline. Each checkout uses its own declared Cargo features.

## Results (2026-09-18)

Header/threading implementation: `52f23c2a`. macOS 26.5.1 arm64, rustc 1.95.0
(59807616e 2026-04-14), release opt-level 3/LTO/one codegen unit. Input is
`hela_gt20peps.mzspeclib.txt.gz`, 266,335,688 bytes, SHA-256
`d58e749ad193fcd97dda63734b18b0fc6edd4e77499f0dc1b734721449d8d636`.

Three passes per mode, in this order:

1. main, branch legacy, serial record, 2 workers, 4 workers;
2. 4 workers, 2 workers, serial record, branch legacy, main;
3. serial record, main, 2 workers, branch legacy, 4 workers.

These measurements predate configurable batches: workers handed off one record
at a time. The preliminary pilot is excluded.

| Mode | Three elapsed times (s) | Median (s) | Correct split |
| --- | --- | --- | --- |
| Serial record | 14.553287, 14.679114, 14.571128 | 14.571128 | Yes |
| Record: 2 workers + reader | 15.156578, 15.155520, 15.210726 | 15.156578 | Yes |
| Record: 4 workers + reader | 8.887028, 8.934134, 9.076343 | 8.934134 | Yes |
| Branch legacy | 32.285228, 32.125165, 33.415145 | 32.285228 | Yes |
| Actual main 8524b24e | 27.975092, 28.468860, 27.876043 | 27.975092 | No |

Four workers give **1.63×** the serial-record throughput and
**3.61×** branch-legacy throughput. Two workers are slightly
slower than serial in this setup; thread dispatch is not free and scaling is not
linear. These timings alone do not identify the limiting stage. The gzip reader
and record framing remain serial; worker scheduling/allocator effects were not
profiled.

Every branch-legacy/record pass matched every integer counter, including all 26
residue bins. Main's combined counters also matched, but its split was incorrect
(all 948,957 spectra assigned to targets).

| Group | Spectra/analytes | Residues | Carbon atoms |
| --- | ---: | ---: | ---: |
| Target | 474,630 | 7,884,870 | 39,529,386 |
| Decoy | 474,327 | 7,882,632 | 39,518,457 |

Amino-acid histogram (letters with zero counts in both groups omitted):

| Amino acid | Target | Decoy |
| --- | ---: | ---: |
| A | 549,765 | 549,585 |
| C | 108,516 | 108,513 |
| D | 435,495 | 435,438 |
| E | 682,911 | 682,563 |
| F | 259,320 | 259,275 |
| G | 431,979 | 431,790 |
| H | 167,556 | 167,544 |
| I | 382,905 | 382,872 |
| K | 599,226 | 599,010 |
| L | 784,440 | 784,266 |
| M | 309,114 | 308,934 |
| N | 293,889 | 293,859 |
| P | 360,057 | 359,979 |
| Q | 396,177 | 396,108 |
| R | 461,193 | 460,839 |
| S | 525,396 | 525,258 |
| T | 383,895 | 383,841 |
| V | 481,950 | 481,896 |
| W | 73,815 | 73,806 |
| Y | 197,271 | 197,256 |

Verification: 158 library unit tests, eight integration tests, three doctests and
three chemistry-example tests passed at that revision (the opt-in snapshot capture remained ignored).
The scoped-record test verifies Send/Sync bounds, shared malformed-header value
isolation and raw/peak buffer reuse after returning from a worker.

## Reusable batch comparison (2026-09-18)

The configurable-batch version keeps the same worker-local reduction. Only buffer
handoffs change: one batch travels to a worker and back, retaining its vector and
record allocations. The CLI defaults to one record; pass `32` explicitly for the
configuration below. No parser code changes were needed.

A six-configuration pilot tried 2/4 workers and batches of 1/8/32. Four workers
with batches of 32 were fastest in that pilot. A separate confirmation used four
workers in batch order 1, 32, 32, 1, 1, 32, with the same direct-gzip input and timer
boundaries as above:

| Records per handoff | Three elapsed times (s) | Median (s) | Reusable record slots |
| --- | --- | ---: | ---: |
| 1 | 8.961180, 8.769287, 9.467465 | 8.961180 | 4 |
| 32 | 6.021132, 6.089533, 6.076543 | 6.076543 | 128 |

Batching gives **1.47× throughput** (32% less elapsed time) in this comparison.
Every target/decoy counter and histogram bin matched the earlier verified results
exactly. The tradeoff is more retained buffers; peak memory was not measured.
This does not establish an optimal batch size for other libraries or workloads.

## Optional allocation instrumentation

Build with `--features allocation-counting` to install a counting wrapper around
`std::alloc::System` in this example only. With the feature disabled, neither the
wrapper nor its atomic counter updates are compiled into the executable. No new
dependency or library allocator is introduced.

```sh
cargo build --release --offline --locked -p mzannotate --example library_chemistry_benchmark --features allocation-counting
target/release/examples/library_chemistry_benchmark legacy 0 ~/fasta/hela_gt20peps.mzspeclib.txt.gz
target/release/examples/library_chemistry_benchmark legacy 4 ~/fasta/hela_gt20peps.mzspeclib.txt.gz 1
target/release/examples/library_chemistry_benchmark legacy 4 ~/fasta/hela_gt20peps.mzspeclib.txt.gz 32
target/release/examples/library_chemistry_benchmark record 0 ~/fasta/hela_gt20peps.mzspeclib.txt.gz
target/release/examples/library_chemistry_benchmark record 4 ~/fasta/hela_gt20peps.mzspeclib.txt.gz 1
target/release/examples/library_chemistry_benchmark record 4 ~/fasta/hela_gt20peps.mzspeclib.txt.gz 32
```

Legacy parallel mode uses the same bounded worker/return-channel schedule as the
record mode, with worker-local counters merged at join. It recycles the outer
`Vec<AnnotatedSpectrum>`; clearing it drops each spectrum's internal allocations.
The reader still eagerly materializes spectra. Both APIs support batches of 1 or
32 with four workers, keeping at most 4 or 128 spectrum slots respectively.

Counters start after shared ontology initialization, before file/header opening,
and finish after worker joins and reader cleanup, before printing:

- `alloc_calls`: successful fresh allocations, including zeroed allocations.
- `free_calls`: deallocations; reallocations are counted separately.
- `realloc_calls`: successful reallocations, including those performed in place.
- `alloc_requested_bytes`: sum of fresh allocation sizes.
- `realloc_requested_bytes`: sum of full requested destination sizes, not bytes copied.
- `growth_bytes`: fresh sizes plus positive realloc size increases; a churn metric.
- `baseline_live_bytes`, `end_live_bytes`, `peak_live_bytes`: live requested Rust
  heap bytes, including allocations retained from initialization. Peak is reset at
  measurement start. Subtract baseline to report additional peak requested bytes.

These are Rust allocator requests, not all native malloc calls. They exclude
allocator metadata, fragmentation, thread stacks and transient storage inside
realloc; peak live requested bytes is not process RSS. Concurrent callbacks make
peak accounting an approximation of the physical instantaneous peak. Counters do
not measure allocator CPU time or bytes actually copied.

**Use a feature-disabled build for speed measurements.** Global atomic counters
add contention and can distort parallel scaling. Run one process per mode; begin/
finish assume no unrelated concurrent allocation traffic. The counter test checks
fresh/zeroed allocations, realloc growth/shrink, free and peak accounting using an
independent allocator instance.

Current example validation: four tests pass without the feature; eight with it.
Legacy/record counters agree across 1/2/4 workers and batches of 1/2/8/32 on a
nine-spectrum fixture, exercising partial batches and reuse. Worker and producer
errors terminate cleanly.

### Allocation results (2026-09-18)

One instrumented release pass per mode on the same 948,957-spectrum gzip file.
All six target/decoy residue histograms and carbon totals match the verified
uninstrumented outputs exactly. Sizes below use decimal GB/MB. Allocation counts
exclude reallocations, shown separately. Peak additional heap is peak requested
live bytes minus the initialization baseline (~142.249 MB), not RSS.

| API | Workers | Batch | Allocations | Reallocations | Growth (GB) | Additional peak (MB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| record | 0 | 1 | 218,354,683 | 13,605,134 | 15.250 | 0.733 |
| record | 4 | 1 | 218,354,748 | 13,605,173 | 15.250 | 0.788 |
| record | 4 | 32 | 218,355,616 | 13,606,785 | 15.251 | 2.344 |
| legacy | 0 | 1 | 498,535,631 | 21,077,344 | 93.165 | 0.813 |
| legacy | 4 | 1 | 498,535,675 | 21,077,344 | 93.165 | 0.993 |
| legacy | 4 | 32 | 498,535,675 | 21,077,344 | 93.165 | 3.310 |

The record workload uses about 230 fresh allocations/spectrum versus 525 for
legacy materialization. Its cumulative growth is about 6.1× lower. This includes
skipping unused peak/fragment decoding as well as record reuse; it does not isolate
the contribution of reuse alone. Batching changes allocation counts very little,
while increasing retained live storage. Batching reduces handoff
frequency; these counters do not attribute CPU time.

Exact byte accounting and frees (same runs):

| API / workers / batch | Frees | Fresh requested bytes | Realloc destination bytes | Growth bytes | Baseline live | End live | Peak live |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| record / 0 / 1 | 218354564 | 13910227386 | 2805223025 | 15250063822 | 142249278 | 142341254 | 142982483 |
| record / 4 / 1 | 218354628 | 13910243705 | 2805268136 | 15250102758 | 142249278 | 142341302 | 143037051 |
| record / 4 / 32 | 218355496 | 13910586713 | 2807055359 | 15251341510 | 142249279 | 142341303 | 144593677 |
| legacy / 0 / 1 | 498535512 | 87567217663 | 11179439595 | 93165001640 | 142249278 | 142341254 | 143062685 |
| legacy / 4 / 1 | 498535555 | 87567225815 | 11179439595 | 93165009792 | 142249278 | 142341302 | 143242552 |
| legacy / 4 / 32 | 498535555 | 87567254583 | 11179439595 | 93165038560 | 142249279 | 142341303 | 145558961 |

Instrumented elapsed times are deliberately omitted from speed comparisons.

### Stage audit before the reuse correction (63f94561)

The aggregate chemistry workload is not an allocation test of record framing
alone. A serial stage probe skipped 100 warm-up spectra, then measured the next
10,000 from the same HeLa gzip input. Each stage used allocator counter deltas;
chemical results were consumed before advancing the record.

| Stage | Fresh allocations / record | Reallocations / record | Growth bytes / record |
| --- | ---: | ---: | ---: |
| Refill record | 0 | 0.0002 | 0.3424 |
| Spectrum/analyte metadata and typed values | 0 | 0.0001 | 0.0768 |
| `analytes()` | 44.0963 | 8.4298 | 6394.8557 |
| Residue/formula counting after cached analytes | 184.5095 | 5.2045 | 9654.6152 |
| Numeric peaks and raw annotation views | 0 | 0.0006 | 1.0368 |

An independent extra probe of ProForma decoding measured 17.0115 allocations and
2.348 reallocations per record; **it is already part of `analytes()`**, not another
stage to add to its cost. The sample differs from the full-library distribution.
Rare reallocations in the view stages reflect buffers encountering larger records.

At that revision, `analytes()` aliased owned legacy `Analyte`: its `owned_groups` path
formats and reparses borrowed metadata, allocates protein/parameter objects and
then discards nested storage on reset. Formula calculation also constructs owned
intermediates. Retaining the outer vectors does **not** establish chemistry reuse.
The PR must not be described as end-to-end allocation-free parsing/chemistry.

A feature-gated regression test now checks zero allocations, frees **and**
reallocations over 100 refills of an already-warmed fixture while traversing all
scope metadata/typed values and numeric/raw peak views, including EOF invalidation.
Test-only thread-local counters exclude other test threads. This proves the fixed
fixture path; it does not cover chemical decoding, structured errors, new maximum
sizes or all possible libraries.

## Reuse correction

Implementation: `33331c38`, with a subsequent guard preserving general-resolver
label ordering and overflow behavior for general formula modifications. Changes are in the
parser/chemistry path, not just in the benchmark:

- `analytes()` no longer builds owned protein/custom metadata through `owned_groups`.
  Those properties remain borrowed scope views until explicit materialization.
- Nested analyte and interpretation slots retain capacity when scopes shrink,
  sections fail or chemistry switches between peptide, formula and unknown.
- `ProFormaScratch` reuses peptide/residue/modification/charge storage for the
  common linear syntax. Molecular-formula text parses into retained element storage.
- `DecodedAnalyte::formulas()` caches borrowed formula alternatives. Unambiguous
  linear formulas accumulate directly into retained buffers instead of allocating
  an intermediate formula collection per residue. The old target-level `.formulas()`
  method remains available and allocating.

The reusable ProForma path is conservative: uppercase linear sequences, numeric
UNIMOD/MOD side-chain modifications, optional nonzero integer proton charges, and
matching placement rules. Other syntax and diagnostic-producing cases use the
full parser. Global isotopes, ambiguous residues/modifications, glycans, cross-links
and labelled/floating-mass formula cases use the general formula resolver when
needed. Those fallbacks, fragment resolution, owned exports and errors can still
allocate. This is not a universal zero-allocation guarantee for all chemistry.

### Full-library allocation counts after correction

One release/instrumented run per configuration; all six complete the original
948,957-spectrum gzip workload and match every target/decoy count and histogram.
Counts include opening, buffers, cache warm-up, worker setup and cleanup. No lookup table of previously seen spectra or ProForma strings is retained; storage is bounded by the reusable
record slots and the largest shapes encountered. Instrumented times are not used
for throughput comparisons.

| API | Workers | Batch | Fresh allocations | Reallocations | Growth bytes | Additional peak requested bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| record | 0 | 1 | 239 | 43 | 743,880 | 740,488 |
| record | 4 | 1 | 418 | 106 | 823,216 | 819,632 |
| record | 4 | 32 | 5,998 | 2,692 | 3,729,360 | 3,721,808 |
| legacy | 0 | 1 | 498,535,631 | 21,077,344 | 93,165,001,640 | 813,407 |
| legacy | 4 | 1 | 498,535,675 | 21,077,344 | 93,165,009,792 | 956,859 |
| legacy | 4 | 32 | 498,535,675 | 21,077,344 | 93,165,038,560 | 6,945,402 |

Serial record access fell from **218,354,683 to 239 fresh allocations**, and from
15,250,063,822 to **743,880 cumulative growth bytes**. Four-worker batch-32 storage
has higher initialization/growth cost because 128 records retain their own buffers.
Peak is requested live Rust heap above the initialized baseline (~142.249 MB), not
RSS; it varies with scheduling. All these counters cover the whole run, not just
a warmed steady-state window.

### Correctness and reuse checks

- Every one of the **948,957 analytes** was independently compared with the original
  ProForma decoder, including reported charge application: full targets, all formula
  alternatives and all formula labels matched.
- All **12 legacy fixture snapshots** match both `63f94561` built with the current
  lockfile and the original characterization outputs. The original snapshot
  command refuses the older lockfile fingerprint; this separate matched-lockfile
  comparison leaves that baseline untouched.
- Warmed allocation tests cover complete metadata/peak/chemistry/formula access,
  alternating long/short scopes, peptide/formula/unknown transitions, interpretation
  reference buffers, and good → bad → good chemistry. Valid warmed decoding performs
  zero allocations, frees or reallocations. Disposal of an owned error is outside
  the post-error decoding measurement.
- Core differential tests cover modifications, charges, warnings, malformed input,
  global isotopes, ambiguous residues, cross-links, complex adducts, formula parsing
  and custom formula labels. Conservative fallbacks preserve original results.
- Full suites passed: mzannotate 158 library tests + 8 integration tests; mzcore
  667 library tests (10 ignored), plus integration/doc tests. The added label-order
  and intermediate-overflow regressions passed in the four targeted reuse tests.
  All eight feature-enabled chemistry-example tests passed, including post-error reuse.

No chemistry-counting policy changed: all analytes and chains are counted, formula
ambiguity remains an error for this consumer, and target/decoy policy is unchanged.

### Uninstrumented timings after correction

Release build without `allocation-counting`, implementation `2c994c08`; subsequent
changes add an intermediate-overflow fallback and move modification handles instead
of cloning them. These timings precede those changes. Same gzip file, timer boundaries
and warm file cache as above.
Three passes, with all integer outputs verified after every run:

1. record serial, legacy 4/batch-32, record 4/batch-32, legacy serial;
2. reverse order;
3. record 4/batch-32, record serial, legacy serial, legacy 4/batch-32.

The third legacy-serial run overlapped a validation build; it was discarded and
repeated without build traffic. Instrumented runs are excluded.

| API | Workers | Batch | Three times (s) | Median (s) |
| --- | ---: | ---: | --- | ---: |
| record | 0 | 1 | 6.600318, 6.544904, 6.527653 | 6.544904 |
| record | 4 | 32 | 5.548407, 5.529183, 5.545259 | 5.545259 |
| legacy | 0 | 1 | 31.289644, 31.725339, 30.900169 | 31.289644 |
| legacy | 4 | 32 | 24.956230, 24.997719, 24.841362 | 24.956230 |

The new four-worker workload is **4.50× faster than batched legacy workers**.
Serial record processing is **4.78× faster than serial legacy**. The new four-worker
path is only 1.18× faster than new serial processing: after reducing chemistry work,
more workers provide a smaller benefit. These timings do not attribute CPU costs
to decompression versus framing. The comparison includes reusable formula calculation
and skipped peak/fragment materialization; it does not isolate text-parser speed.
