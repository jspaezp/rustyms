# Temporary compatibility snapshots

From the repository root, before changing the parser:

```sh
python3 scripts/mzspeclib_snapshots.py baseline
python3 scripts/mzspeclib_snapshots.py check
```

The harness discovers every **tracked** `.mzspeclib.txt`/`.mzspeclib` input, including
the fuzz corpus, case-insensitively. Generated `.out` files are excluded. It captures
the header, explicit materialized fields, serialized chemical structures and peaks,
structured errors, iteration through failures, repeated EOF, `read_next` at EOF,
and indexed reads by index/key/name/literal scan number. Identical indexed results
reference sequential events to avoid duplicating peak arrays.

Unordered maps are sorted. Analyte `params` and `proteins` are also canonicalized
as multisets: the current reader populates them by randomized HashMap traversal,
and repeat captures confirmed their order differs within unchanged runs. Their
order is therefore not covered; multiplicity and values remain covered. Every
other array retains its order. Fixture-relative paths avoid checkout-specific diagnostics.
This is not a write/reparse round-trip or a Debug snapshot. Existing Serde formats
are used for supporting types; dependency changes require separate review. These
fixtures do not exhaust malformed-input, I/O-error, or nonfinite-number behavior;
targeted permanent regression tests cover intentional fixes and new contracts.

Artifacts live in gitignored `.mzspeclib-snapshots/`. The manifest records the base
commit, uncommitted parser-source diff, fixture hashes, capture-code hash, compiler,
and Cargo.lock hash; the lockfile is retained alongside snapshots. Baseline creation
refuses overwrites. Checks reject fixture, lockfile or capture-code changes, keep
candidate outputs and readable diffs, and never bless mismatches automatically.

Run a check immediately after capture to detect nondeterminism before refactoring.
If existing behavior is nondeterministic, document and narrowly address the cause;
do not sort meaningful arrays just to pass. Intentional bug fixes require focused
tests and a separately reviewed behavior change, not blanket snapshot regeneration.
Keep the baseline until the PR's compatibility review is complete, then delete the
generated directory. No bulk snapshot artifacts belong in Git history.
