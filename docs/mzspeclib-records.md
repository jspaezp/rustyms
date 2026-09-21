# Lazy mzSpecLib records

The record API preserves supplied metadata and peak text before projecting anything
into `AnnotatedSpectrum`. Numeric peaks, analyte chemistry, interpretations, and
annotation reports are independently requested batches. Repeated access borrows
cached results; `.iter()` performs no decoding or collection allocation.

```rust
use mzannotate::mzspeclib::record::MzSpecLibLibrary;
use mzcore::ontology::STATIC_ONTOLOGIES;

let text = b"<mzSpecLib>\n<Spectrum=1>\n<Peaks>\n100\t1\t?\n";
let mut library = MzSpecLibLibrary::open(&text[..], None, &STATIC_ONTOLOGIES)?;
let mut reader = library.reader();
let mut record = reader.empty_record();
while reader.read_into(&mut record)? {
    for attribute in record.attributes()?.iter() {
        println!("{} = {}", attribute.name(), attribute.raw_value());
    }
    for peak in record.peaks()?.iter() {
        println!("{} {} {:?}", peak.mz(), peak.intensity(), peak.annotation_field());
    }
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

Use `reader.records()` to retain independently owned records. They borrow the
library's header/ontology context, not the reader handle. The library must outlive
them. `read_into` instead reuses raw text, indexes, flat numeric/annotation arrays,
and cache capacities, including across parse errors. EOF or failure clears the
current record; check its return value before accessing it again.

Views contain references and ranges. They can outlive temporary collection views,
but Rust prevents replacing a record while its views are still used. Records are
`Send` but not `Sync`: each worker owns its record, whose lazy caches use safe
`OnceCell`/`RefCell` storage. The borrowed header/context is immutable and `Sync`. Raw peak text is still
read and buffered during metadata-only traversal. Decoded analyte/interpretation
slots retain nested storage across refills, smaller scopes, target-kind changes and
failed chemistry. Common linear chemistry and formula accumulation reuse storage
after warm-up. Complex chemistry, fragment resolution, owned exports, diagnostics
and newly encountered larger shapes can still allocate; this is not a universal
zero-allocation guarantee.

## Framing without producer-side parsing

For producer/worker pipelines, use `empty_frame()` and `read_frame_into()` instead
of `empty_record()` and `read_into()`. A `SpectrumFrame` owns the reusable record
storage; it does not copy the text when exposing a parsed record.

```rust
let mut reader = library.reader();
let mut frame = reader.empty_frame();
while reader.read_frame_into(&mut frame)? {
    // A pipeline sends the frame to its worker at this point.
    // No IDs, attributes or peak fields have been parsed by the producer.
    let record = frame.record()?;
    for attribute in record.attributes()?.iter() {
        println!("{} = {}", attribute.name(), attribute.raw_value());
    }
    // Return the frame to the producer after these borrowed views expire.
}
```

`read_frame_into` resets raw framing state, finds `<Spectrum=` line boundaries,
and buffers the complete text with source coordinates. Decoded-cache cleanup is
deferred to the next `frame.record()` call on the consumer. IO and invalid UTF-8 remain
reader errors. `frame.record()` parses structural metadata once on the calling
thread, retains its index, and returns a borrowed `SpectrumRecord`. Chemistry,
metadata values, inheritance, peaks and annotations retain their existing lazy
accessors. Structural errors are cached and leave raw text/source coordinates
available for diagnostics until refill. EOF or IO failure invalidates the frame.

The frame is `Send`, not `Sync`. Borrowed record views prevent refill, so consumers
cannot accidentally use stale indexes. Text, indexes and nested chemistry buffers
all survive a worker round trip. `read_into` and `records()` remain convenience
APIs that prepare structural metadata before returning; their existing error
behavior is unchanged. Use frames when that work belongs on a worker.

## Consumer-owned early rejection

Run the executable dogfood example from the repository root:

```sh
cargo run -p mzannotate --example library_loading -- \
  data/fetal_brain_tiny_consensus_td.mzSpecLib.txt
cargo run -p mzannotate --example library_loading -- \
  data/fetal_brain_tiny_consensus_td.mzSpecLib.txt --peaks
cargo run -p mzannotate --example library_loading -- \
  data/fetal_brain_tiny_consensus_td.mzSpecLib.txt --annotations
```

The example inspects effective spectrum-origin occurrences and rejects the shipped
fixture's shuffle-and-reposition decoys before requesting peaks or chemistry. Its
explicit policy recognizes only `MS:1003195`; production ancestry classification
belongs in the consumer. Unknown/absent origin is retained as unclassified, not
converted to a target assertion. Predicted and decoy assertions can coexist.

Previously this required either fully parsing a spectrum and then rejecting it,
or a raw preflight pass that duplicated attribute-set resolution. Metadata access
now succeeds even when unrequested numeric peaks or ProForma would fail.

## Evidence and diagnostics

| Need | API |
| --- | --- |
| Distinct identities | `key()`, declared-name occurrences, `source_position()` |
| Raw/effective attributes | `local_attributes()`, `attributes()` |
| Inherited provenance | occurrence `origin()` includes set/claim/source group |
| Defaults despite overrides | `header().attribute_sets()` |
| Undefined claims | `attribute_set_claims()`, `resolution_diagnostics()` |
| Typed term identity | `occurrence.value()` → `ValueView::Term`, namespace-aware `Curie` |
| Original text | `raw_value()`; unsupported accessions retain name/raw spelling |
| Units and custom pairing | `occurrence.group()`, `custom_values()`; inspect ambiguity |
| All analyte metadata without chemistry | `analyte_scopes()` |
| Interpreted chemistry | `analytes()`; original scope evidence remains available |
| Borrowed calculated formulas | `analyte.formulas()`; cached alternatives, no implicit first choice |
| Explicit field presence | `SuppliedField::Absent` versus `Present("")` / `Present("?")` |
| Peak provenance | `source_row()`, `source_position()`, aggregation and extra columns |
| Every annotation alternative | `resolved_annotations()?.iter()`, each row's `iter()` |
| Malformed fields with zero alternatives | report `diagnostics()` |

Header metadata values are classified once during `open`, with errors stored per
occurrence rather than making malformed values fail the whole header. Record-local
metadata values remain lazy: one cached batch per raw record buffer, with outcomes
stored per occurrence. An unrelated malformed value does not hide a
valid typed origin. Lexical views never copy text. Numeric conversions such as
`to_f64()` are fallible; no missing coordinate is replaced with zero. Quantity
selection, unit normalization, ontology ancestry, competition groups and fatality
of unknown evidence remain consumer policy.

`Ok(annotation_report)` means the report was constructed, not that every annotation
resolved. Always inspect `diagnostics()`: it includes whole-field errors even when
that row has no parsed alternatives. Invalid analyte chemistry is a soft dependency
failure: raw fields and recoverable annotation syntax remain in the report, with
blocked-resolution diagnostics. Numeric peak parsing is a hard prerequisite.
Reported deviation and peak m/z stay separate from optional resolved fragment data.

Interpretation-member scopes remain raw, explicitly unsupported evidence; strict
`validate()`/materialization rejects them. Unknown attribute-set claims remain
inspectable diagnostics. `materialize()` creates an owned, strict projection;
`into_annotated()` moves decoded analytes/interpretations/fragments into that output.
Neither is an alternative error-policy mode for the legacy reader.

Source IDs and effective group IDs are opaque process-local ownership tokens, not
persistent IDs or ordinals. Record generations change on reuse. A group ID from a
different record, scope or previous generation cannot match accidentally. Persist
the supplied spectrum key plus your file identity instead.

## Indexed access

For `BufRead + Seek` inputs beginning at byte zero, call `build_index()` and
`read_by_index_into`, `read_by_key_into`, `read_by_name_into`, or
`read_by_scan_number_into`. Indexing scans declarations/literal spectrum metadata;
it does not decode chemistry or peaks. Building the index preserves reader position.
Lookups work after EOF and reset lookahead/caches before loading. Name and scan
lookups intentionally index literal spectrum attributes, not inherited names.

## Legacy compatibility

`MzSpecLibTextParser` retains its existing implementation and public behavior.
Its source-order chemistry/error handling is kept deliberately: forcing new section
caches would change which competing error appears first and would read farther
before errors. The new API is additive. It applies consistent inheritance to
analytes/interpretations; those corrections do not silently change legacy output.

The temporary snapshot harness covers all 12 tracked fixture libraries, including
fuzz copies. See [snapshot instructions](../mzannotate/tests/mzspeclib/SNAPSHOTS.md).
The only canonicalized arrays are analyte proteins/parameters already produced in
nondeterministic HashMap order. Snapshots are local ignored artifacts, not bundled
test fixtures. No universal compatibility claim is made for inputs absent from the
corpus; targeted tests cover new contracts and recovery/reuse behavior.

## Aggregate workload benchmark

On the supplied HeLa library (948,957 spectra), summing intensities and averaging
reported fragment/selected-ion m/z by target/decoy took 6.217 s with reusable records
versus 22.981 s with the legacy reader (median of three passes; decompression excluded).
Counts and aggregates matched. See [benchmark source, methodology and results](mzspeclib-benchmark.md).

Direct gzip loading takes 8.956 s versus 25.241 s (same workload,
three-pass medians). Actual main (`8524b24e`) takes 22.166 s but incorrectly
places every spectrum in the target group; see the benchmark's separate main results.

## Parallel chemistry with reusable records

[library_chemistry_benchmark.rs](../mzannotate/examples/library_chemistry_benchmark.rs)
uses scoped threads and a reusable batch per worker (default: one record). The
reader moves a batch of raw `SpectrumFrame`s to its worker; after indexing/decoding/counting, the worker
returns ownership for refill. Counters remain worker-local and merge at `join`. Raw text, numeric/metadata cache storage and cache-vector capacity move
with the record. No record/raw-buffer clones or Arc-wrapped library context are
needed. Chemical objects can still allocate during decoding.

The library must outlive the scoped workers. Individual records are not shared
concurrently. Gzip decompression and record framing stay on the reader thread;
structural indexing, metadata resolution, analyte decoding and carbon counting run on workers. The
bounded channels cap reusable records at worker count × batch size. Worker errors terminate processing instead of losing records
or emitting partial counts as successful results.

See [chemistry benchmark methodology and results](mzspeclib-chemistry-benchmark.md)
for serial/parallel amino-acid and carbon counting on the supplied HeLa library.

## Chemistry reuse and owned export

`DecodedAnalyte` contains an ID and decoded target, with private reusable chemistry
and formula buffers. It is no longer an alias for the owned legacy `Analyte`.
Protein/custom properties stay in `analyte_scopes().attributes()`; chemistry access
does not stringify and reparse those properties. Likewise, decoded interpretations
contain IDs, probability and reusable analyte-reference arrays; their other
properties stay in interpretation scope views. `materialize()`/`into_annotated()`
explicitly create the owned legacy protein/parameter/attribute structures.

```rust
for analyte in record.analytes()?.iter() {
    // Same cached slice on repeated access; all alternatives remain available.
    for formula in analyte.formulas().iter() {
        inspect_elements(analyte.id, formula.elements());
    }
}
```

`analyte.target.formulas()` remains the old allocating chemistry trait method.
Use **`analyte.formulas()`** for the record's reusable, cached calculation. Unknown
targets return an empty formula slice, reported molecular formulas return one, and
ambiguous chemistry retains every alternative. Declared formula text remains
available through scope attributes, separately from calculated formulas.

The reusable ProForma decoder handles uppercase linear residues with numeric
`UNIMOD:`/`MOD:` side-chain modifications and optional nonzero integer proton
charge. It uses the existing ontology and placement predicates. Named ions,
terminal modifications, global isotope syntax, glycans, cross-links, ambiguity,
complex adducts and diagnostic-producing inputs go through the full parser;
those paths retain semantics but may allocate. No unsupported syntax is silently
simplified. Formula accumulation avoids intermediates for unambiguous linear
chemistry with formula/database modifications; general cases use the existing
resolver. Fragment-annotation chemistry has not acquired this zero-allocation
common path.

The warmed allocation tests include modified peptides, metadata/peak views,
interpretation references, changing scope counts, peptide/formula/unknown target
transitions, and successful chemistry after an error. They check allocation,
reallocation **and free** calls; owned diagnostic disposal is outside the
post-error decoding measurement. See the [chemistry benchmark](mzspeclib-chemistry-benchmark.md)
for full-file measurements and limits.
