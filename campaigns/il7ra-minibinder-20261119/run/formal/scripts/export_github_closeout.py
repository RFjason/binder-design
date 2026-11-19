"""Collect immutable completed-campaign evidence on a Slurm CPU node."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import zipfile

assert os.environ.get("SLURM_JOB_ID"), "Run through Slurm"
ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT.with_name("il7ra_minibinder_20261110")
OUT = ROOT / "publication/il7ra-closeout-20261119"
assert not OUT.exists(), "Use a fresh output directory"
GIT = OUT / "git"
EVIDENCE = OUT / "evidence"


def load(path):
    return json.loads(path.read_text())


def dump(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def copy(src, dest):
    assert src.is_file(), str(src)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)


summary = load(ROOT / "reports/completed.json")
assert summary["total_scored"] == 2855 and summary["passing_total"] == 445
rows = [json.loads(line) for line in (ROOT / "reports/scored_candidates.jsonl").read_text().splitlines()]
passing = sorted((r for r in rows if r["pilot_gate"]), key=lambda r: (-r["ipSAE_min"], r["id"]))
assert len(passing) == 445
shortlist = load(ROOT / "reports/combined_shortlist.json")
assert len(shortlist) == 20
raw_by_id = {r["id"]: (r, origin) for origin in (PILOT, ROOT)
             for r in load(origin / "reports/candidates.json")}

for name in ("completed.json", "scored_candidates.csv", "scored_candidates.jsonl",
             "diversity.json", "pilot_audit.json", "combined_shortlist.json"):
    copy(ROOT / "reports" / name, GIT / "results" / name)
copy(ROOT / "reports/shortlist.csv", GIT / "top20/candidates.csv")
copy(ROOT / "reports/shortlist.fasta", GIT / "top20/candidates.fasta")
for label, origin in (("formal", ROOT), ("pilot", PILOT)):
    for folder in ("inputs", "scripts", "config"):
        for src in sorted((origin / folder).glob("*")):
            if src.is_file() and src.suffix in (".json", ".fasta", ".a3m", ".pdb", ".cif", ".py", ".slurm", ".tsv"):
                copy(src, GIT / "run" / label / folder / src.name)
    for name in ("scoring_protocol.json",):
        if (origin / name).is_file():
            copy(origin / name, GIT / "run" / label / name)
    for name in ("chosen_patches.json", "patch_geometry.json", "preflight.json", "final_preflight.json",
                 "collection_wave_0.json", "collection_wave_1.json", "collection_wave_2.json"):
        if (origin / "reports" / name).is_file():
            copy(origin / "reports" / name, GIT / "run" / label / "reports" / name)
    for name in ("controller.json", "closeout.json", "jobs.json"):
        if (origin / "state" / name).is_file():
            copy(origin / "state" / name, GIT / "run" / label / "state" / name)
    for name in ("3DI3.pdb", "anthropic_results_v1_cached.csv"):
        if (origin / "sources" / name).is_file():
            copy(origin / "sources" / name, GIT / "run" / label / "sources" / name)
    for i in range(1, 13):
        src = origin / "predictions" / f"control_{i:02d}" / "summary.json"
        if src.is_file():
            copy(src, GIT / "run" / label / "controls" / f"control_{i:02d}.json")

provenance = []
for index, row in enumerate(passing, 1):
    cid = row["id"]
    raw, origin = raw_by_id[cid]
    assert raw["sequence"] == row["sequence"]
    record = {"id": cid, "sequence": row["sequence"], "origin": "formal" if origin == ROOT else "pilot",
              "sources": raw["sources"], "prediction_files": {}}
    for kind in ("complex", "monomer"):
        srcdir = origin / "predictions" / (cid + "_" + kind)
        assert load(srcdir / "summary.json")["ok"]
        for src in sorted(srcdir.iterdir()):
            if src.is_file():
                copy(src, EVIDENCE / cid / kind / src.name)
        for name in ("model.pdb", "summary.json"):
            copy(srcdir / name, GIT / "passing/structures" / cid / kind / name)
        record["prediction_files"][kind] = f"evidence/{cid}/{kind}"
    source_files = []
    for i, source in enumerate(raw["sources"]):
        src = origin / source["file"]
        digest = hashlib.sha256(src.read_bytes()).hexdigest()
        assert digest == source["source_sha256"], str(src)
        dest = EVIDENCE / cid / "generation" / f"{i:02d}_{src.name}"
        copy(src, dest)
        source_files.append({"original_file": source["file"], "archive_path": str(dest.relative_to(OUT)), "sha256": digest})
    record["source_files"] = source_files
    provenance.append(record)
    if index % 50 == 0:
        print(f"COLLECTED {index}/{len(passing)}", flush=True)

for row in shortlist:
    srcdir = ROOT / "predictions" / (row["id"] + "_historical")
    assert load(srcdir / "summary.json")["ok"]
    for src in sorted(srcdir.iterdir()):
        if src.is_file():
            copy(src, EVIDENCE / row["id"] / "historical" / src.name)
    for name in ("model.pdb", "summary.json"):
        copy(srcdir / name, GIT / "top20/historical" / row["id"] / name)

dump(GIT / "provenance/passing_sources.json", provenance)
archive = OUT / "il7ra-passing-445-full-evidence.zip"
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as handle:
    for path in sorted(EVIDENCE.rglob("*")):
        if path.is_file():
            handle.write(path, path.relative_to(OUT))
with zipfile.ZipFile(archive) as handle:
    assert handle.testzip() is None
    entries = len(handle.namelist())
dump(GIT / "provenance/archive_manifest.json", {"release_tag": "il7ra-closeout-20261119", "parts": [{
    "filename": archive.name, "bytes": archive.stat().st_size,
    "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(), "entries": entries,
    "candidate_count": 445, "historical_rechecks": 20,
    "contents": "Full complex/monomer prediction directories, original generation structures, and Top20 historical predictions."}]})
bundle = OUT / "il7ra-github-data.tar.gz"
import tarfile
with tarfile.open(bundle, "w:gz", compresslevel=6) as handle:
    handle.add(GIT, arcname="git")
print(json.dumps({"status": "EXPORTED", "candidates": len(passing), "top20": len(shortlist),
                  "evidence_bytes": archive.stat().st_size, "git_bundle_bytes": bundle.stat().st_size,
                  "slurm_job_id": os.environ["SLURM_JOB_ID"]}), flush=True)
