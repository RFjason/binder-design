"""Select crystal-supported patches and make paired pilot inputs."""

from itertools import combinations
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np
from biotite import structure as struc
from biotite.structure.io.pdb import PDBFile

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
ROOT = Path(__file__).resolve().parents[1]
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
              "ARNDCQEGHILKMFPSTWYV"))


def structure(path):
    atoms = PDBFile.read(path).get_structure(model=1, extra_fields=["b_factor"])
    return atoms[atoms.element != "H"]


def protein(atoms, chain):
    return atoms[(atoms.chain_id == chain) & np.isin(atoms.res_name, list(AA))]


def ca_by_position(atoms):
    return {int(a.res_id) + 20: a.coord for a in atoms[atoms.atom_name == "CA"]}


def fit(source, target):
    ids = sorted(source.keys() & target.keys())
    a = np.array([source[i] for i in ids])
    b = np.array([target[i] for i in ids])
    ac, bc = a.mean(axis=0), b.mean(axis=0)
    u, _, vt = np.linalg.svd((a - ac).T @ (b - bc))
    correction = np.diag([1, 1, np.linalg.det(u @ vt)])
    rotation = u @ correction @ vt
    translation = bc - ac @ rotation
    rmsd = np.sqrt(np.mean(np.sum((a @ rotation + translation - b) ** 2, axis=1)))
    return rotation, translation, float(rmsd)


def distance(a, b):
    return float(np.min(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1))) if len(b) else 999.0


def residue_surface(atoms, sasa=None):
    if sasa is None:
        sasa = struc.sasa(atoms)
    return {int(i) + 20: float(np.nansum(sasa[(atoms.res_id == i) &
                   ~np.isin(atoms.atom_name, ["N", "C", "O", "CA"])]))
            for i in np.unique(atoms.res_id)}


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


complex_ = structure(ROOT / "sources/3DI3.pdb")
bound = protein(complex_, "B")
apo = protein(structure(ROOT / "sources/3UP1.pdb"), "A")
preview = "--crystal-preview" in sys.argv
canonical = (ROOT / "inputs/IL7RA_E21_D239.fasta").read_text().splitlines()[1]
model_dir = ROOT / "predictions/target_E21_D239"
if preview:
    model = apo
else:
    assert json.loads((model_dir / "summary.json").read_text())["ok"]
    model = protein(structure(model_dir / "model.pdb"), "A")
    model_ca = model[model.atom_name == "CA"]
    assert "".join(AA[a] for a in model_ca.res_name) == canonical
    assert list(model_ca.res_id) == list(range(1, 220))
    assert all(set(model.atom_name[model.res_id == i]) >= {"N", "CA", "C", "O"} for i in model_ca.res_id)
model_ca = model[model.atom_name == "CA"]
confidence_scale = 100 if max(model_ca.b_factor) <= 1.1 else 1
confidence = {int(a.res_id) + 20: float(a.b_factor * confidence_scale) for a in model_ca}
il7 = protein(complex_, "A")
glycans = complex_[np.isin(complex_.res_name, ["NAG", "BMA", "MAN"])]
assert len(glycans) >= 60
model_ca_by = ca_by_position(model)
bound_ca_by, apo_ca_by = ca_by_position(bound), ca_by_position(apo)
rotation, translation, bound_rmsd = fit(bound_ca_by, model_ca_by)
_, _, apo_rmsd = fit(apo_ca_by, model_ca_by)
assert max(bound_rmsd, apo_rmsd) <= 3.5, (bound_rmsd, apo_rmsd)
glycan_in_model = glycans.coord @ rotation + translation
bound_sasa, apo_sasa = residue_surface(bound), residue_surface(apo)
model_atom_sasa = struc.sasa(model)
model_sasa = residue_surface(model, model_atom_sasa)
sequons = [i + 21 for i in range(len(canonical) - 2)
           if canonical[i] == "N" and canonical[i+1] != "P" and canonical[i+2] in "ST"]
asn = model.coord[np.isin(model.res_id + 20, sequons) & (model.atom_name == "ND2")]
ends = model.coord[(model.res_id + 20 < 37) | (model.res_id + 20 > 229)]
rows = []
for i in sorted(bound_ca_by.keys() & apo_ca_by.keys()):
    mb = model[model.res_id == i - 20]
    xb = bound[bound.res_id == i - 20]
    side_indices = np.flatnonzero((model.res_id == i - 20) &
                                  ~np.isin(model.atom_name, ["N", "C", "O", "CA"]) &
                                  (model_atom_sasa >= 4))
    if not len(side_indices) or i == 66 or canonical[i-21] in "CGP":
        continue
    exposed_atoms = sorted(side_indices, key=lambda j: -model_atom_sasa[j])[:2]
    d_il7 = distance(xb.coord, il7.coord)
    clearance = min(distance(xb.coord, glycans.coord), distance(mb.coord, glycan_in_model))
    rows.append({"position": i, "residue": canonical[i-21], "il7_distance": d_il7,
                 "glycan_distance": clearance, "sequon_nd2_distance": distance(mb.coord, asn),
                 "termini_distance": distance(mb.coord, ends),
                 "surface": min(bound_sasa[i], apo_sasa[i], model_sasa[i]),
                 "plddt": None if preview else confidence[i],
                 "atom_names": [str(model.atom_name[j]) for j in exposed_atoms]})


def patch(mode):
    pool = [r for r in rows if r["surface"] >= 24 and (preview or r["plddt"] >= 75) and
            r["glycan_distance"] >= 8 and r["sequon_nd2_distance"] >= 8 and
            r["termini_distance"] >= 10 and
            (r["il7_distance"] <= 5 if mode == "il7_face" else r["il7_distance"] >= 13)]
    triples = []
    for group in combinations(pool, 3):
        ca = [model_ca_by[r["position"]] for r in group]
        span = max(np.linalg.norm(a-b) for a, b in combinations(ca, 2))
        if span <= 14:
            quality = sum(r["surface"] + 1.5 * min(r["glycan_distance"], 20) for r in group)
            triples.append((quality - 2*span, group))
    if not triples:
        raise RuntimeError(f"No exposed contiguous triple for {mode}. Inspect reports/patch_geometry.json.")
    return list(max(triples, key=lambda pair: pair[0])[1])


reports = ROOT / "reports"
suffix = "_crystal_preview" if preview else ""
dump(reports / f"patch_geometry{suffix}.json", {"bound_model_ca_rmsd": bound_rmsd, "apo_model_ca_rmsd": apo_rmsd,
                                             "sequons": sequons, "positions": rows})
patches = {mode: patch(mode) for mode in ("il7_face", "distal_face")}
dump(reports / f"chosen_patches{suffix}.json", patches)
print(json.dumps({"bound_rmsd": bound_rmsd, "apo_rmsd": apo_rmsd, "patches": patches}), flush=True)
if preview:
    raise SystemExit(0)

inputs = ROOT / "inputs"
shutil.copyfile(model_dir / "model.pdb", inputs / "receptor.pdb")
shutil.copyfile(model_dir / "model.cif", inputs / "receptor.cif")
batches = []
for mode, hotspots in patches.items():
    numbers = [r["position"] - 20 for r in hotspots]
    for low, high in [(55, 75), (76, 95)]:
        for method, repeats in [("rfd4", 1), ("bc2", 6)]:
            for replica in range(repeats):
                name = f"{method}_{mode}_{low}_{replica:02d}"
                seed = 711000 + len(batches)
                out = ROOT / "generation" / name
                config = inputs / (name + ".json")
                batches.append({"id": name, "method": method, "patch": mode, "min_length": low,
                                "max_length": high, "seed": seed, "replicates": 12 if method == "rfd4" else 1})
                if method == "rfd4":
                    dump(config, {name: {"input": str(inputs / "receptor.cif"),
                         "contig": f"{low}-{high}, /0, A1-219",
                         "conditions": {"C_CRDGRP": {"1": [{"select": "sel('A/*/*/[N,CA,C,O]')"}]},
                           "C_SEQ": {"true": [{"select": "sel('A')"}]},
                           "C_HOT": {"true": [{"select": f"sel('A/*/{number}/[{','.join(row['atom_names'])}]')"}
                                             for number, row in zip(numbers, hotspots)]}}}})
                else:
                    dump(config, {
                         "targets": [{"name": "IL7RA", "target_path": str(inputs / "receptor.pdb"),
                                      "chains": "A", "hotspots": ",".join(map(str, numbers))}],
                         "binder_lengths": [low, high], "number_of_final_designs": 1, "max_trajectories": 1,
                         "auto_multi_gpu": False, "design_workers": 1, "campaign_seed": seed,
                         "min_plddt_refine": 0.65, "min_iptm_anneal": 0.5, "min_iptm_harden": 0.5,
                         "min_plddt_mutate": 0.65, "min_iptm_mutate": 0.5, "min_iptm_final": 0.6,
                         "min_monomer_plddt_final": 0.7, "project_folder": str(out), "aa_bias": {"C": 0},
                         "save_failed_refolds": True, "save_failed_trajectories": True,
                         "archive_trajectories": False})
dump(inputs / "design_batches.json", batches)
assert len(batches) == 28
