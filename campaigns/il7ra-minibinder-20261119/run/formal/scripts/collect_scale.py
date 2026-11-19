"""Register unique candidates serially and score only unfinished predictions."""
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
ROOT = Path(__file__).resolve().parents[1]
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
              "ARNDCQEGHILKMFPSTWYV"))
wave = int(os.environ["SCALE_WAVE"])


def read(path, default=None):
    return json.loads(path.read_text()) if path.is_file() else default


def dump(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


config = read(ROOT / "config/campaign.json")
batches = read(ROOT / "inputs/design_batches.json")
target = "".join(line.strip() for line in (ROOT / "inputs/IL7RA_E21_D239.fasta").read_text().splitlines()
                 if not line.startswith(">"))
baseline = set(read(ROOT / "state/baseline_sequences.json"))
with (ROOT / "sources/anthropic_results_v1_cached.csv").open(newline="") as handle:
    public = {row["sequence"] for row in csv.DictReader(handle)
              if row["target"].lower().replace("-", "") == "il7ra"}
assert len(public) >= 12
records = {r["sequence"]: r for r in read(ROOT / "reports/candidates.json", [])}
counts, rejected, duplicates, over_quota = {}, [], [], []
sources = {}
for batch in batches:
    if batch["wave"] != wave:
        continue
    folder = ROOT / "generation" / batch["id"]
    if batch["method"] == "rfd4":
        files = sorted((folder / "inference_samples").rglob("*_generated.cif"))
    else:
        files = sorted(p for p in folder.rglob("*.cif") if "Complexes" in p.parts)
        files += sorted((folder / "1_Trajectories").rglob("trajectory.cif"))
    counts[batch["id"]] = {"source_files": len(files),
        "generation_done": (ROOT / "reports" / ("done_" + batch["id"] + ".json")).is_file()}
    sources.setdefault(batch["patch"], []).extend((batch, path) for path in files)

# Interleave patch groups so the final quota cannot be filled by one patch's files alone.
pending = []
while any(sources.values()):
    for files in sources.values():
        if files:
            pending.append(files.pop(0))
for batch, path in pending:
    try:
        atoms = get_structure(CIFFile.read(path), model=1)
        ca = atoms[(atoms.atom_name == "CA") & np.isin(atoms.res_name, list(AA))]
        sequences = {str(c): "".join(AA[x] for x in ca.res_name[ca.chain_id == c])
                     for c in np.unique(ca.chain_id)}
        target_ids = [cid for cid, sequence in sequences.items() if sequence == target]
        assert len(target_ids) == 1 and len(sequences) == 2, "Target sequence or chain count mismatch"
        target_id = target_ids[0]
        binder_id = next(cid for cid in sequences if cid != target_id)
        sequence = sequences[binder_id]
        assert batch["min_length"] <= len(sequence) <= batch["max_length"]
        assert np.isfinite(ca.coord).all()
        assert sequence not in public, "Exact public sequence"
        if sequence in baseline:
            duplicates.append(str(path.relative_to(ROOT)))
            continue
        if sequence not in records and len(records) >= config["new_quota"]:
            over_quota.append(str(path.relative_to(ROOT)))
            continue
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        parent = batch["id"] + ":" + (path.stem if batch["method"] == "rfd4" else "trajectory_1")
        source = {**batch, "parent": parent, "file": str(path.relative_to(ROOT)),
                  "target_chain": target_id, "binder_chain": binder_id,
                  "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        record = records.setdefault(sequence, {"id": "IL7RA_s5k_" + digest[:16],
                                               "sequence": sequence, "sources": []})
        if not any(s["file"] == source["file"] for s in record["sources"]):
            record["sources"].append(source)
    except Exception as exc:
        rejected.append({"file": str(path.relative_to(ROOT)), "error": repr(exc)})

candidates = list(records.values())
dump(ROOT / "reports/candidates.json", candidates)
cases = []
for row in candidates:
    binder = {"id": "B", "sequence": row["sequence"]}
    for suffix, chains in (("complex", [{"id": "A", "sequence": target,
                                         "msa": "inputs/IL7RA_E21_D239.a3m"}, binder]),
                           ("monomer", [binder])):
        case_id = row["id"] + "_" + suffix
        if read(ROOT / "predictions" / case_id / "summary.json", {}).get("ok"):
            continue
        cases.append({"id": case_id, "chains": chains})
for shard in range(config["score_shards"]):
    dump(ROOT / "inputs" / f"wave_{wave}_score_{shard}.json", cases[shard::config["score_shards"]])
report = {"wave": wave, "registered_new": len(candidates), "baseline": config["baseline_count"],
          "unfinished_prediction_cases": len(cases), "counts": counts, "rejected": rejected,
          "baseline_duplicates": duplicates, "over_quota_files": len(over_quota)}
dump(ROOT / "reports" / f"collection_wave_{wave}.json", report)
print(json.dumps({k: report[k] for k in ("wave", "registered_new", "unfinished_prediction_cases")}), flush=True)
