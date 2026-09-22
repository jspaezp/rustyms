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

## Reproduce the decompressed-file comparison

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

For the decompressed-file comparison, build time, process startup, ontology
initialization, decompression and result printing are outside the reported timers. Opening/reading the decompressed file,
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
excluded from every parse timing above. Direct streaming gzip timings are measured separately below; they cannot be
inferred by simply adding these preparation and parser timings.

## Direct gzip loading, including current main

The original "legacy" comparison above used the existing reader on the feature
branch, which includes the attribute-set fixes from its parent PR. It did **not**
measure `main`. The follow-up checks freshly fetched `origin/main` commit
`8524b24ebfe62fa308532e61a8a12c4b3ae522d5` in a detached worktree, with no parser
source changes. Only the benchmark example and its `flate2` dev dependency are
added to that checkout. The legacy aggregation loop is copied unchanged; imports
and the match arm referencing the unavailable record API are removed.

Main has no committed Cargo.lock. Its build is seeded with the feature checkout's
lockfile, then resolved offline. This removes unused dependencies but introduces
no new registry package versions. Each checkout retains its own Cargo features.
Both use the same compiler and release profile.

All three paths read the original compressed file directly, using `flate2`'s
`MultiGzDecoder`, a 256 KiB compressed-input buffer and a 256 KiB decoded-text
buffer. The timer includes gzip decoding, file reads, parsing and aggregation.
Build/process startup, ontology initialization and output printing remain excluded.
Runs are sequential with warm file cache.

```sh
cargo build --release --offline --locked -p mzannotate --example library_benchmark
target/release/examples/library_benchmark record ~/fasta/hela_gt20peps.mzspeclib.txt.gz
target/release/examples/library_benchmark legacy ~/fasta/hela_gt20peps.mzspeclib.txt.gz
```

Main also completed one decompressed-file pass in **19.770023 s**, but its grouping
was wrong: every spectrum was assigned to targets. Empty-group means print `NaN`.
This is the result exposed by main's metadata API under the same classifier, not
a raw-input reclassification or a consumer workaround for missing inheritance.

To prepare the same isolated main executable from the repository root:

```sh
benchmark_main_dir=$(mktemp -d)
git worktree add --detach "$benchmark_main_dir" 8524b24ebfe62fa308532e61a8a12c4b3ae522d5
cp Cargo.lock "$benchmark_main_dir/Cargo.lock"
uv run --no-project python - "$benchmark_main_dir" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1])
source = Path('mzannotate/examples/library_benchmark.rs').read_text()
source = source.replace(
    'use mzannotate::mzspeclib::{\n    MzSpecLibTextParser,\n    record::{MzSpecLibLibrary, ValueView},\n};',
    'use mzannotate::mzspeclib::MzSpecLibTextParser;',
).replace('use mzcv::curie;\n', '')
start = source.index('        "record" => {')
end = source.index('        _ => return Err("mode must be legacy or record".into()),', start)
source = source[:start] + source[end:]
example = root / 'mzannotate/examples/library_benchmark_main.rs'
example.parent.mkdir(exist_ok=True)
example.write_text(source)
manifest = root / 'mzannotate/Cargo.toml'
manifest.write_text(manifest.read_text() + '\n[dev-dependencies]\nflate2 = { workspace = true }\n')
PY
CARGO_TARGET_DIR="$PWD/target" cargo build \
  --manifest-path "$benchmark_main_dir/Cargo.toml" \
  --release --offline -p mzannotate --example library_benchmark_main
target/release/examples/library_benchmark_main legacy ~/fasta/hela_gt20peps.mzspeclib.txt.gz
git worktree remove --force "$benchmark_main_dir"
```

### Direct gzip results

Order: record, branch legacy, main, main, branch legacy, record, record,
branch legacy, main. Three runs per path:

| Implementation | Times (s) | Median (s) | Correct target/decoy split |
| --- | --- | --- | --- |
| Record API | 8.995650, 8.955765, 8.944596 | 8.955765 | Yes |
| Branch legacy | 25.240764, 25.291547, 25.177315 | 25.240764 | Yes |
| Main 8524b24e | 22.175208, 22.131116, 22.166280 | 22.166280 | No |

Record API is **2.82× faster** than branch legacy for direct
gzip loading with matching results. Main's timing is not an equivalent-correctness
comparison: all three main runs produced this incorrect split:

```text
target spectra=948957 peaks=27287703 intensity_sum=4830421.133137887344 fragment_mz_mean=812.217455227002 precursor_mz_mean=693.135861399064
decoy spectra=0 peaks=0 intensity_sum=0.000000000000 fragment_mz_mean=NaN precursor_mz_mean=NaN
```
