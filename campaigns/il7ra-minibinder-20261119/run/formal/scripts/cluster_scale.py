"""Cluster passing candidates by TM-align fold and receptor contact footprint."""
import json
import os
from pathlib import Path
import sys

import numpy as np
from biotite.structure.io.pdb import PDBFile

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".vendor"))
from tmtools import tm_align


def load(path):
    return json.loads(path.read_text())


audit = load(ROOT / "reports/pilot_audit.json")
config = load(ROOT / "config/campaign.json")
baseline = load(Path(config["baseline_audit"]))
rows = [(r, ROOT) for r in audit["rows"] if r["pilot_gate"]]
rows += [(r, Path(config["baseline_root"])) for r in baseline["rows"] if r["pilot_gate"]]
rows.sort(key=lambda pair: -pair[0]["ipSAE_min"])
fold_leaders, interface_leaders = [], []
members, fold_members = [], {}
for row, origin in rows:
    atoms = PDBFile.read(origin / "predictions" / (row["id"] + "_monomer") / "model.pdb").get_structure(model=1)
    coords = np.asarray(atoms.coord[(atoms.chain_id == "B") & (atoms.atom_name == "CA")], dtype=np.float64)
    assert len(coords) == len(row["sequence"])
    fold = None
    for index, (leader_coords, leader_sequence) in enumerate(fold_leaders):
        score = tm_align(coords, leader_coords, row["sequence"], leader_sequence)
        if min(score.tm_norm_chain1, score.tm_norm_chain2) >= 0.65:
            fold = index
            break
    if fold is None:
        fold = len(fold_leaders)
        fold_leaders.append((coords, row["sequence"]))
    if "target_contacts" in row:
        contacts = set(row["target_contacts"])
    else:
        complex_ = PDBFile.read(origin / "predictions" / (row["id"] + "_complex") / "model.pdb").get_structure(model=1)
        receptor = complex_[(complex_.chain_id == "A") & (complex_.element != "H")]
        binder = complex_[(complex_.chain_id == "B") & (complex_.element != "H")]
        distances = np.linalg.norm(receptor.coord[:, None] - binder.coord[None, :], axis=-1)
        contacts = {int(r) + 20 for r in receptor.res_id[distances.min(axis=1) <= 5]}
    interface = None
    for index, leader in enumerate(interface_leaders):
        if len(contacts & leader) / max(1, len(contacts | leader)) >= 0.6:
            interface = index
            break
    if interface is None:
        interface = len(interface_leaders)
        interface_leaders.append(contacts)
    sources = [{**s, "wave": s.get("wave", -1),
                "patch": "known" if origin != ROOT else s["patch"]} for s in row["passing_sources"]]
    members.append({"id": row["id"], "fold_cluster": fold, "interface_cluster": interface,
                    "ipSAE_min": row["ipSAE_min"], "sources": sources, "origin": str(origin)})
    fold_members.setdefault(str(fold), []).append(row["id"])

group_stats = {}
for group in ("known", "adjacent", "other"):
    parents, joint_first = {}, {}
    for member in members:
        matching = [s for s in member["sources"] if s["patch"].split("_")[0] == group]
        for source in matching:
            parents[source["parent"]] = source["wave"]
            key = (member["fold_cluster"], member["interface_cluster"])
            joint_first[key] = min(joint_first.get(key, source["wave"]), source["wave"])
    group_stats[group] = {"passing_parents": len(parents), "parent_waves": list(parents.values()),
                          "joint_cluster_first_waves": list(joint_first.values()),
                          "joint_clusters": len(joint_first)}
report = {"wave": int(os.environ["SCALE_WAVE"]), "passing_sequences": len(rows), "fold_clusters": len(fold_leaders),
          "interface_clusters": len(interface_leaders), "members": members, "groups": group_stats,
          "fold_min_both_normalizations": 0.65, "interface_jaccard_threshold": 0.6}
path = ROOT / "reports/diversity.json"
temporary = path.with_suffix(".json.tmp")
temporary.write_text(json.dumps(report, indent=2) + "\n")
temporary.replace(path)
print(json.dumps({k: report[k] for k in ("passing_sequences", "fold_clusters", "interface_clusters", "groups")}),
      flush=True)
