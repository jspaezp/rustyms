#!/usr/bin/env python3
"""Capture/check temporary mzSpecLib characterization snapshots (stdlib only)."""
import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / ".mzspeclib-snapshots"


def run(*args):
    return subprocess.check_output(args, cwd=ROOT).decode().strip()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixtures():
    paths = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    return sorted(p for p in paths if p.lower().endswith((".mzspeclib.txt", ".mzspeclib")))


def manifest(paths):
    return {
        "schema": 1,
        "commit": run("git", "rev-parse", "HEAD"),
        "rustc": run("rustc", "--version"),
        "fixtures": {p: digest(ROOT / p) for p in paths},
        "cargo_lock_sha256": digest(ROOT / "Cargo.lock"),
        "parser_diff": run("git", "diff", "HEAD", "--", "mzannotate/src", "mzcore/src", "mzcv/src"),
        "capture_sha256": digest(ROOT / "mzannotate/tests/mzspeclib/snapshots.rs"),
    }


def capture(output, paths):
    output.mkdir(parents=True)
    config = output / "config.json"
    config.write_text(json.dumps({"root": str(ROOT), "output": str(output), "fixtures": paths}))
    env = dict(os.environ, MZSPECLIB_SNAPSHOT_CONFIG=str(config))
    subprocess.run([
        "cargo", "test", "--offline", "--locked", "-p", "mzannotate", "--test", "mzspeclib",
        "snapshots::capture_fixture_snapshots", "--", "--ignored", "--exact", "--nocapture",
    ], cwd=ROOT, env=env, check=True, timeout=600)
    config.unlink()
    (output / "manifest.json").write_text(json.dumps(manifest(paths), indent=2) + "\n")
    # Preserve the dependency resolution for reproducing this baseline later.
    (output / "Cargo.lock").write_bytes((ROOT / "Cargo.lock").read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["baseline", "check"])
    parser.add_argument("--baseline", type=Path, default=ARTIFACTS / "baseline")
    args = parser.parse_args()
    baseline = args.baseline.resolve()
    paths = fixtures()
    if not paths:
        parser.error("no tracked mzSpecLib fixtures found")
    ARTIFACTS.mkdir(exist_ok=True)
    if args.mode == "baseline":
        if baseline.exists():
            parser.error(f"refusing to overwrite baseline: {baseline}")
        with tempfile.TemporaryDirectory(prefix="capture-", dir=ARTIFACTS) as staging:
            output = Path(staging) / "output"
            capture(output, paths)
            baseline.parent.mkdir(parents=True, exist_ok=True)
            output.rename(baseline)
        print(f"Captured {len(paths)} fixtures in {baseline}")
        return
    expected = json.loads((baseline / "manifest.json").read_text())
    current = manifest(paths)
    for field in ("schema", "fixtures", "cargo_lock_sha256", "capture_sha256"):
        if current[field] != expected[field]:
            parser.error(f"{field} changed; comparison would not isolate parser behavior")
    # Retain failed comparisons for inspection; never modify the baseline.
    output = Path(tempfile.mkdtemp(prefix="check-", dir=ARTIFACTS)) / "output"
    capture(output, paths)
    changed = []
    for fixture in paths:
        name = fixture + ".json"
        before = (baseline / name).read_text()
        after = (output / name).read_text()
        if before != after:
            changed.append(fixture)
            diff = output / (name + ".diff")
            diff.write_text("".join(difflib.unified_diff(
                before.splitlines(True), after.splitlines(True),
                fromfile=f"baseline/{fixture}", tofile=f"current/{fixture}",
            )))
            print(f"Changed: {fixture}; diff: {diff}")
    if changed:
        raise SystemExit(f"{len(changed)}/{len(paths)} fixture snapshots changed")
    print(f"All {len(paths)} fixture snapshots match baseline {expected['commit']}; output: {output}")


if __name__ == "__main__":
    main()
