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
read and buffered during metadata-only traversal. Molecular objects, owned exports,
and error construction can allocate; this is not a zero-allocation chemistry parser.

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
reader moves a loaded batch to its worker; after decoding/counting, the worker
returns ownership for refill. Counters remain worker-local and merge at `join`. Raw text, numeric/metadata cache storage and cache-vector capacity move
with the record. No record/raw-buffer clones or Arc-wrapped library context are
needed. Chemical objects can still allocate during decoding.

The library must outlive the scoped workers. Individual records are not shared
concurrently. Gzip decompression and record framing stay on the reader thread;
metadata resolution, analyte decoding and carbon counting run on workers. The
bounded channels cap reusable records at worker count × batch size. Worker errors terminate processing instead of losing records
or emitting partial counts as successful results.

See [chemistry benchmark methodology and results](mzspeclib-chemistry-benchmark.md)
for serial/parallel amino-acid and carbon counting on the supplied HeLa library.
