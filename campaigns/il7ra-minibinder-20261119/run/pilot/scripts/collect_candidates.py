"""Extract both tools' de novo designs and build balanced scoring shards."""
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
inputs = ROOT / "inputs"
receptor = {"id": "A", "sequence": (inputs / "IL7RA_E21_D239.fasta").read_text().splitlines()[1],
            "msa": "inputs/IL7RA_E21_D239.a3m"}
controls = json.loads((inputs / "anthropic_controls_reference_only.json").read_text())
with (ROOT / "sources/anthropic_results_v1_cached.csv").open(newline="") as handle:
    public = [row for row in csv.DictReader(handle) if row["target"].lower().replace("-", "") == "il7ra"]
assert len(public) >= len(controls), "Public reference table target mismatch"
public_sequences = {row["sequence"] for row in public} | {c["sequence"] for c in controls}
batches = json.loads((inputs / "design_batches.json").read_text())
records = {}
counts, rejected = {}, []
for batch in batches:
    source = ROOT / "generation" / batch["id"]
    if batch["method"] == "rfd4":
        files = sorted((source / "inference_samples").rglob("*_generated.cif"))
    else:
        files = sorted(p for p in source.rglob("*.cif") if "Complexes" in p.parts)
        files += sorted((source / "1_Trajectories").rglob("trajectory.cif"))
    counts[batch["id"]] = {"files": len(files), "accepted_sources": 0,
                            "job_complete": (ROOT / "reports" / ("done_" + batch["id"] + ".json")).exists()}
    for path in files:
        try:
            atoms = get_structure(CIFFile.read(path), model=1)
            ca = atoms[(atoms.atom_name == "CA") & np.isin(atoms.res_name, list(AA))]
            sequences = {str(c): "".join(AA[x] for x in ca.res_name[ca.chain_id == c])
                         for c in np.unique(ca.chain_id)}
            target_ids = [c for c in sequences if sequences[c] == receptor["sequence"]]
            assert len(target_ids) == 1 and len(sequences) == 2, "Target sequence or chain count mismatch"
            target_chain = target_ids[0]
            binder_chain = next(c for c in sequences if c != target_chain)
            seq = sequences[binder_chain]
            assert batch["min_length"] <= len(seq) <= batch["max_length"], "Unexpected binder length"
            assert seq not in public_sequences, "Exact public/control sequence"
            assert np.isfinite(ca.coord).all(), "Nonfinite structure"
            parent = batch["id"] + ":" + (path.stem if batch["method"] == "rfd4" else "trajectory_1")
            entry = {**batch, "parent": parent, "file": str(path.relative_to(ROOT)),
                     "target_chain": target_chain, "binder_chain": binder_chain,
                     "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            digest = hashlib.sha256(seq.encode()).hexdigest()
            record = records.setdefault(seq, {"id": "IL7RA_pilot_" + digest[:16], "sequence": seq, "sources": []})
            record["sources"].append(entry)
            counts[batch["id"]]["accepted_sources"] += 1
        except Exception as exc:
            rejected.append({"file": str(path.relative_to(ROOT)), "reason": repr(exc)})


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


candidates = list(records.values())
dump(ROOT / "reports/collection.json", {"counts": counts, "rejected": rejected,
                                        "public_reference_count": len(public), "unique": len(candidates)})
dump(ROOT / "reports/candidates.json", candidates)
assert candidates, "No valid candidates found. See collection report and generation logs."
for i in range(16):
    cases = []
    for row in candidates[i::16]:
        binder = {"id": "B", "sequence": row["sequence"]}
        cases.extend([{"id": row["id"] + "_complex", "chains": [receptor, binder]},
                      {"id": row["id"] + "_monomer", "chains": [binder]}])
    dump(inputs / f"pilot_score_{i}.json", cases)
print(f"COLLECTED {len(candidates)} unique candidates from {sum(r['files'] for r in counts.values())} files",
      flush=True)
