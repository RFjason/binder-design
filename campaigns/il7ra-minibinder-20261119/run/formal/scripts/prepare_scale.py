"""Prepare the 5,000-candidate campaign from the completed IL7RA pilot."""
import ast
from itertools import combinations
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
from biotite.structure.io.pdb import PDBFile
from esm.utils.msa import MSA

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT.with_name("il7ra_minibinder_20261110")


def read(path):
    return json.loads(path.read_text())


def dump(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


for name in ("inputs", "sources", "reports", "predictions", "config", "state", "cache", "logs"):
    (ROOT / name).mkdir(exist_ok=True)
assert not (ROOT / "config/campaign.json").exists(), "Existing campaign must not be overwritten"
for name in ("IL7RA_E21_D239.fasta", "IL7RA_E21_D239.a3m", "IL7RA_competition_M1_G236.fasta",
             "IL7RA_competition_M1_G236.a3m", "anthropic_controls_reference_only.json",
             "receptor.cif", "receptor.pdb"):
    shutil.copyfile(PILOT / "inputs" / name, ROOT / "inputs" / name)
for name in ("3DI3.pdb", "anthropic_results_v1_cached.csv"):
    shutil.copyfile(PILOT / "sources" / name, ROOT / "sources" / name)
shutil.copyfile(PILOT / "reports/patch_geometry.json", ROOT / "reports/patch_geometry.json")
shutil.copyfile(PILOT / "scoring_protocol.json", ROOT / "scoring_protocol.json")
for index in range(1, 13):
    name = f"control_{index:02d}"
    folder = ROOT / "predictions" / name
    folder.mkdir(exist_ok=True)
    shutil.copyfile(PILOT / "predictions" / name / "summary.json", folder / "summary.json")

protocol = read(ROOT / "scoring_protocol.json")
sequence = "".join(line.strip() for line in (ROOT / protocol["target"]["fasta"]).read_text().splitlines()
                   if not line.startswith(">"))
assert len(sequence) == 219 and hashlib.sha256(sequence.encode()).hexdigest() == protocol["target"]["sequence_sha256"]
msa_path = ROOT / protocol["target"]["msa"]["file"]
assert hashlib.sha256(msa_path.read_bytes()).hexdigest() == protocol["target"]["msa"]["sha256"]
msa = MSA.from_a3m(msa_path, remove_insertions=True)
assert msa.depth == 1033 and msa.sequences[0] == sequence
for path in (ROOT / "scripts").glob("*.py"):
    ast.parse(path.read_text(), filename=str(path))

atoms = PDBFile.read(ROOT / "inputs/receptor.pdb").get_structure(model=1)
ca = atoms[(atoms.chain_id == "A") & (atoms.atom_name == "CA")]
assert list(ca.res_id) == list(range(1, 220))
points = {int(a.res_id) + 20: a.coord for a in ca}
geometry = read(ROOT / "reports/patch_geometry.json")
safe = [r for r in geometry["positions"] if r["surface"] >= 24 and r["plddt"] >= 75
        and r["glycan_distance"] >= 8 and r["sequon_nd2_distance"] >= 8 and r["termini_distance"] >= 10]
known = read(PILOT / "reports/chosen_patches.json")["il7_face"]


def center(group):
    return np.array([points[r["position"]] for r in group]).mean(axis=0)


def span(group):
    return max(np.linalg.norm(points[a["position"]] - points[b["position"]])
               for a, b in combinations(group, 2))


def quality(group):
    return sum(r["surface"] + 1.5 * min(r["glycan_distance"], 20) for r in group) - 2 * span(group)


patches = {"known": known}
known_positions = {r["position"] for r in known}
known_center = center(known)
adjacent_options = []
for group in combinations([r for r in safe if r["il7_distance"] <= 12], 3):
    shift = np.linalg.norm(center(group) - known_center)
    if span(group) <= 14 and 5 <= shift <= 16 and len(known_positions & {r["position"] for r in group}) <= 1:
        adjacent_options.append(group)
for name in ("adjacent_1", "adjacent_2"):
    eligible = [g for g in adjacent_options if all(
        np.linalg.norm(center(g) - center(patches[p])) >= 5
        and len({r["position"] for r in g} & {r["position"] for r in patches[p]}) <= 1
        for p in patches if p.startswith("adjacent"))]
    assert eligible, f"No geometrically eligible patch for {name}"
    patches[name] = list(max(eligible, key=quality))

other_options = [g for g in combinations([r for r in safe if r["il7_distance"] >= 13], 2)
                 if span(g) <= 12 and np.linalg.norm(center(g) - known_center) >= 18]
for name in ("other_1", "other_2"):
    eligible = [g for g in other_options if all(
        np.linalg.norm(center(g) - center(patches[p])) >= 10
        for p in patches if p.startswith("other"))]
    assert eligible, f"No geometrically eligible patch for {name}"
    patches[name] = list(max(eligible, key=quality))
dump(ROOT / "reports/chosen_patches.json", patches)

prior = read(PILOT / "reports/candidates.json")
assert len(prior) == 78
dump(ROOT / "state/baseline_sequences.json", [r["sequence"] for r in prior])
config = {"target_count": 5000, "baseline_count": 78, "new_quota": 4922,
          "baseline_root": str(PILOT), "baseline_audit": str(PILOT / "reports/pilot_audit.json"),
          "primary_waves": 5, "score_shards": 60, "generation_concurrency": {"rfd4": 4, "bc2": 2},
          "score_concurrency": 2, "rfd4_partition": "a800", "bc2_partition": "a100",
          "score_partition": "a100", "branch_weights": {"known": 0.5, "adjacent": 0.3, "other": 0.2}}


def make_batch(method, patch, low, wave, replica):
    batch_id = f"{method}_{patch}_{low}_w{wave:02d}_{replica:03d}"
    high = 75 if low == 55 else 95
    index = len(batches)
    item = {"id": batch_id, "method": method, "patch": patch, "wave": wave,
            "min_length": low, "max_length": high, "seed": 730000 + index,
            "replicates": 20 if method == "rfd4" else 1,
            "min_hotspots": 1 if patch.startswith("other") else 2}
    hot = patches[patch]
    numbers = [r["position"] - 20 for r in hot]
    if method == "rfd4":
        spec = {batch_id: {"input": str(ROOT / "inputs/receptor.cif"),
            "contig": f"{low}-{high}, /0, A1-219", "conditions": {
                "C_CRDGRP": {"1": [{"select": "sel('A/*/*/[N,CA,C,O]')"}]},
                "C_SEQ": {"true": [{"select": "sel('A')"}]},
                "C_HOT": {"true": [{"select": f"sel('A/*/{n}/[{','.join(r['atom_names'])}]')"}
                                   for n, r in zip(numbers, hot)]}}}}
    else:
        spec = read(PILOT / "inputs/bc2_il7_face_76_00.json")
        spec["targets"][0].update(target_path=str(ROOT / "inputs/receptor.pdb"),
                                  hotspots=",".join(map(str, numbers)))
        spec.update(binder_lengths=[low, high], campaign_seed=item["seed"],
                    project_folder=str(ROOT / "generation" / batch_id))
    dump(ROOT / "inputs" / (batch_id + ".json"), spec)
    batches.append(item)


batches = []
for wave in range(5):
    # Each wave contains 800 RFD4 samples and 60 independent BC2 trajectories.
    for patch, count in (("known", 20), ("adjacent_1", 6), ("adjacent_2", 6), ("other_1", 4), ("other_2", 4)):
        for replica in range(count):
            make_batch("rfd4", patch, 55 if replica % 2 == 0 else 76, wave, replica)
    for patch, count in (("known", 30), ("adjacent_1", 9), ("adjacent_2", 9), ("other_1", 6), ("other_2", 6)):
        for replica in range(count):
            # BC2's proven length gets most trajectories, with a short exploratory fraction.
            low = 55 if patch != "known" and replica % 3 == 0 else 76
            make_batch("bc2", patch, low, wave, replica)
assert sum(b["replicates"] for b in batches if b["method"] == "rfd4") == 4000
assert sum(b["method"] == "bc2" for b in batches) == 300
assert len({b["seed"] for b in batches}) == 500
dump(ROOT / "inputs/design_batches.json", batches)
dump(ROOT / "config/campaign.json", config)
dump(ROOT / "reports/preflight.json", {"ok": True, "patches": {
    name: [r["residue"] + str(r["position"]) for r in group] for name, group in patches.items()},
    "msa_rows": msa.depth, "batches": len(batches), "job": os.environ["SLURM_JOB_ID"]})
print(json.dumps(read(ROOT / "reports/preflight.json")), flush=True)
