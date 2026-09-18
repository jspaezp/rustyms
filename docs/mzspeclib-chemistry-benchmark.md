# Chemistry counting and scoped workers

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

The worker implementation keeps exactly one reusable record per worker. A single
producer decompresses gzip and frames records, then moves each loaded record to a
worker. The worker resolves metadata, decodes analytes, counts chemistry and returns
the record for refill. Worker-local integer counters are merged after completion.
Channels are bounded; the library outlives scoped workers. The input reader stays
on the producer thread and does not need to implement Send. No raw buffers or decoded
analytes are cloned for dispatch. Record caches still use `OnceCell`/`RefCell`; only
the shared header is immutable. Chemical decoding/formula calculation still allocate.
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
results. Example tests cover the modified count and a failing worker without
stranding other workers. The package's record test checks header errors remain
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
import sys
root = Path(sys.argv[1])
s = Path('mzannotate/examples/library_chemistry_benchmark.rs').read_text().split('#[cfg(test)]')[0]
s = s.replace('    record::{MzSpecLibLibrary, SpectrumRecord, ValueView},\n', '').replace('use mzcv::curie;\n', '')
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

The preliminary pilot is excluded.

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
two chemistry-example tests passed (the opt-in snapshot capture remained ignored).
The scoped-record test verifies Send/Sync bounds, shared malformed-header value
isolation and raw/peak buffer reuse after returning from a worker.
