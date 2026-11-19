"""Audit scored pilot candidates and prepare old-construct sensitivity cases."""
import csv
import json
import os
from pathlib import Path
import re

import numpy as np
from biotite import structure as struc
from biotite.sequence import ProteinSequence
from biotite.sequence.align import SubstitutionMatrix, align_optimal
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "predictions"
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
              "ARNDCQEGHILKMFPSTWYV"))


def load(path):
    return json.loads(path.read_text())


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def pdb(path):
    atoms = PDBFile.read(path).get_structure(model=1)
    return atoms[atoms.element != "H"]


def chain(atoms, chain_id):
    return atoms[(atoms.chain_id == chain_id) & np.isin(atoms.res_name, list(AA))]


def ca(atoms):
    return atoms[atoms.atom_name == "CA"]


def pairs(a, b):
    return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)


def contacts(target, binder):
    return set((int(t) + 21, int(b) + 1)
               for t, b in zip(*np.where(pairs(ca(target).coord, ca(binder).coord) <= 8)))


patches = load(ROOT / "reports/chosen_patches.json")
sequons = load(ROOT / "reports/patch_geometry.json")["sequons"]
ref = pdb(ROOT / "sources/3DI3.pdb")
ref_target = ca(chain(ref, "B"))
sugars = ref[np.isin(ref.res_name, ["NAG", "BMA", "MAN"])]
assert len(sugars)
with (ROOT / "sources/anthropic_results_v1_cached.csv").open(newline="") as handle:
    public = [(row["full_name"], row["sequence"]) for row in csv.DictReader(handle)
              if row["target"].lower().replace("-", "") == "il7ra"]
assert len(public) >= 12
matrix = SubstitutionMatrix.std_protein_matrix()
reference_sequences = [(name, ProteinSequence(seq)) for name, seq in public]


def public_hits(sequence):
    query = ProteinSequence(sequence)
    hits = []
    for name, ref_sequence in reference_sequences:
        result = align_optimal(query, ref_sequence, matrix, gap_penalty=(-10, -1),
                               local=True, max_number=1)[0]
        trace = result.trace
        paired = trace[np.all(trace >= 0, axis=1)]
        coverage = len(paired) / len(query)
        identity = sum(query.code[i] == ref_sequence.code[j] for i, j in paired) / max(1, len(trace))
        duplicate = bool(len(paired) / len(ref_sequence) > 0.90 and identity == 1)
        if duplicate or (coverage >= 0.80 and identity > 0.70):
            hits.append({"reference": name, "coverage": coverage, "identity": float(identity),
                         "duplicate": duplicate})
    return hits


def source_geometry(target, binder, source):
    atoms = get_structure(CIFFile.read(ROOT / source["file"]), model=1)
    source_t = chain(atoms, source["target_chain"])
    source_b = chain(atoms, source["binder_chain"])
    _, transform = struc.superimpose(ca(target).coord, ca(source_t).coord)
    aligned_b = transform.apply(ca(source_b).coord)
    pose_rmsd = float(struc.rmsd(ca(binder).coord, aligned_b))
    source_contacts = contacts(source_t, source_b)
    agreement = len(source_contacts & contacts(target, binder)) / max(1, len(source_contacts))
    hot = [r["position"] for r in patches[source["patch"]]]
    t_hot = target[np.isin(target.res_id + 20, hot)]
    d_hot = pairs(t_hot.coord, binder.coord)
    hotspot_hits = sorted({int(x) + 20 for x in t_hot.res_id[d_hot.min(axis=1) <= 5]})
    return {"parent": source["parent"], "patch": source["patch"], "method": source["method"],
            "pose_ca_rmsd_A": pose_rmsd, "design_contact_recall": agreement,
            "hotspots_contacted": hotspot_hits}


rows, failures = [], []
for row in load(ROOT / "reports/candidates.json"):
    case = row["id"]
    try:
        complex_dir, monomer_dir = PRED / (case + "_complex"), PRED / (case + "_monomer")
        score, mono_score = load(complex_dir / "summary.json"), load(monomer_dir / "summary.json")
        complex_ = pdb(complex_dir / "model.pdb")
        target, binder = chain(complex_, "A"), chain(complex_, "B")
        monomer = chain(pdb(monomer_dir / "model.pdb"), "B")
        assert len(ca(target)) == 219 and len(ca(binder)) == len(ca(monomer)) == len(row["sequence"])
        assert list(ca(target).res_id) == list(range(1, 220))
        assert "".join(AA[x] for x in ca(binder).res_name) == row["sequence"]
        _, to_model = struc.superimpose(ca(target).coord[ref_target.res_id - 1], ref_target.coord)
        d_glycan = float(pairs(to_model.apply(sugars.coord), binder.coord).min())
        nd2 = target.coord[np.isin(target.res_id + 20, sequons) & (target.atom_name == "ND2")]
        d_sequon = float(pairs(nd2, binder.coord).min())
        d_end = float(pairs(target.coord[(target.res_id < 17) | (target.res_id > 209)], binder.coord).min())
        clash_pairs = int(np.sum(pairs(target.coord, binder.coord) < 2.2))
        mono_aligned, _ = struc.superimpose(ca(binder).coord, ca(monomer).coord)
        monomer_rmsd = float(struc.rmsd(ca(binder).coord, mono_aligned))
        geometry = []
        for source in row["sources"]:
            try:
                geometry.append(source_geometry(target, binder, source))
            except Exception as exc:
                failures.append({"id": case, "parent": source["parent"], "error": repr(exc)})
        passing = [g for g in geometry if g["pose_ca_rmsd_A"] <= 5 and g["design_contact_recall"] >= 0.4
                   and len(g["hotspots_contacted"]) >= 2]
        cysteines = binder.coord[binder.atom_name == "SG"]
        if len(cysteines):
            distances = pairs(cysteines, cysteines)
            np.fill_diagonal(distances, np.inf)
            paired = bool(np.all(np.sum(distances < 3.0, axis=1) == 1))
        else:
            paired = row["sequence"].count("C") == 0
        preliminary = (score["ipsae_min"] >= 0.6 and mono_score["chain_plddt"]["B"] >= 75
                       and monomer_rmsd <= 4 and d_glycan >= 3.2 and d_sequon >= 4.5
                       and d_end >= 3.2 and clash_pairs == 0 and paired and bool(passing))
        hits = public_hits(row["sequence"]) if preliminary else None
        rows.append({"id": case, "sequence": row["sequence"], "ipSAE_min": score["ipsae_min"],
                     "ipSAE_directions": score["ipsae_directions"], "ipTM": score["ipTM"],
                     "monomer_pLDDT": mono_score["chain_plddt"]["B"], "monomer_rmsd_A": monomer_rmsd,
                     "crystal_glycan_min_A": d_glycan, "sequon_ND2_min_A": d_sequon,
                     "unresolved_end_min_A": d_end, "interchain_clash_pairs": clash_pairs,
                     "cysteine_count": row["sequence"].count("C"), "paired_cysteines": paired,
                     "binder_sequon_count": len(re.findall(r"(?=N[^P][ST])", row["sequence"])),
                     "source_geometry": geometry, "passing_sources": passing,
                     "initial_public_reference_hits": hits, "pilot_gate": preliminary and not hits})
    except Exception as exc:
        failures.append({"id": case, "error": repr(exc)})

rows.sort(key=lambda r: (-r["ipSAE_min"], -r["monomer_pLDDT"], r["id"]))
controls = load(ROOT / "inputs/anthropic_controls_reference_only.json")
control_rows = []
for i, reference in enumerate(controls, 1):
    path = PRED / f"control_{i:02d}/summary.json"
    result = load(path) if path.is_file() else None
    control_rows.append({"name": reference["full_name"], "binder_final": reference["binder_final"],
                         "ipSAE_min": result["ipsae_min"] if result else None})
shortlist, seen_parents = [], set()
for row in rows:
    parent = next((g for g in row["passing_sources"] if g["parent"] not in seen_parents), None)
    if row["pilot_gate"] and parent and len(shortlist) < 20:
        shortlist.append({**row, "selected_source": parent})
        seen_parents.update(g["parent"] for g in row["passing_sources"])
dump(ROOT / "reports/pilot_audit.json", {"scored": len(rows), "failures": failures, "rows": rows,
     "controls": control_rows, "independent_shortlist": shortlist,
     "public_reference_count": len(public), "scope": "Pilot screening. Full glycans and patent novelty pending."})
old = ROOT / "inputs/IL7RA_competition_M1_G236"
old_sequence = "".join(line.strip() for line in old.with_suffix(".fasta").read_text().splitlines()
                       if not line.startswith(">"))
assert len(old_sequence) == 236
old_receptor = {"id": "A", "sequence": old_sequence,
                "msa": "inputs/IL7RA_competition_M1_G236.a3m"}
for i in range(4):
    dump(ROOT / "inputs" / f"historical_score_{i}.json",
         [{"id": row["id"] + "_historical", "historical_recheck": True,
           "chains": [old_receptor, {"id": "B", "sequence": row["sequence"]}]}
          for row in shortlist[i::4]])
print(f"PILOT_AUDIT scored={len(rows)} independent_shortlist={len(shortlist)} failures={len(failures)}", flush=True)
assert rows, "No successfully scored candidates to audit"
