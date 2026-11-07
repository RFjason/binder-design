"""Compute target-anchored historical pose clusters or audit mature-core predictions."""
import argparse
import hashlib
from itertools import combinations
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np
from biotite.structure.io.pdb import PDBFile

ROOT = Path(__file__).resolve().parents[1]
DOCKQ = ROOT.parent / "anthropic_pdl1_pipeline_20261006/envs/scoring/bin/DockQ"
AA = "ACDEFGHIKLMNPQRSTVWY"
THREE = dict(zip("ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(), AA))


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def pdb(path):
    return PDBFile.read(path).get_structure(model=1)


def chain(atoms, cid, sequence):
    atoms = atoms[atoms.chain_id == cid]
    ca = atoms[atoms.atom_name == "CA"]
    assert "".join(THREE[x] for x in ca.res_name) == sequence, f"Unexpected chain sequence: {cid}"
    assert len(set(ca.res_id)) == len(sequence)
    assert np.isfinite(atoms.coord).all()
    return atoms, ca


def frame_match(row, historical=False):
    frame = list(row["framework_core"])
    seq = row["source_sequence"] if historical else row["sequence"]
    if historical and row["binder_changed"]:
        frame[-1] = frame[-1][:-1]
    match = re.fullmatch(f"([{AA}]{{5,30}})".join(map(re.escape, frame)), seq)
    assert match
    return match


def heavy(atoms):
    return atoms[~np.isin(np.char.upper(atoms.element), ("H", "D"))]


def distances(a, b):
    return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)


def aligned_rmsd(a, b):
    # Row-vector Kabsch fit of target B onto A; apply to binder B without a binder fit.
    x, y = a["anchor"], b["anchor"]
    yc, xc = y.mean(axis=0), x.mean(axis=0)
    u, _, vt = np.linalg.svd((y-yc).T @ (x-xc))
    r = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
    delta = (b["frame"]-yc) @ r + xc - a["frame"]
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def cluster(rows, canon):
    anchors = {"crd2": (84, 125), "crd3": (127, 166)}
    poses = {}
    for row in rows:
        cid = row["source_design_id"]
        atoms = pdb(row["source_model"])
        t, tc = chain(atoms, "A", canon[:211])
        b, bc = chain(atoms, "B", row["source_sequence"])
        t = heavy(t[np.isin(t.res_id, tc.res_id[29:])])
        contact = distances(t.coord, heavy(b).coord).min(axis=1) <= 4.5
        indices = {resid: pos for pos, resid in enumerate(tc.res_id, start=1)}
        contacts = {indices[resid] for resid in t.res_id[contact]}
        assert contacts, cid
        site = row["stratum"].split("_")[0]
        lo, hi = anchors[site]
        frame = np.ones(len(bc), dtype=bool)
        for g in (1, 2, 3):
            match = frame_match(row, historical=True)
            frame[match.start(g):match.end(g)] = False
        # The shared scaffold excludes the final, newly restored S128.
        if not row["binder_changed"]:
            frame[-1] = False
        assert frame.sum() == 83
        poses[cid] = dict(anchor=tc.coord[lo-1:hi].astype(float), frame=bc.coord[frame].astype(float),
                          contacts=contacts, site=site, method=row["method"], parent_id=row["parent_id"],
                          ipsae=row["v1_ipsae_min"], cif_sha256=digest(row["source_complex"]),
                          pdb_sha256=digest(row["source_model"]))
    pairs = {}
    for a, b in combinations(sorted(poses), 2):
        x, y = poses[a], poses[b]
        if x["site"] == y["site"]:
            pairs[a, b] = (aligned_rmsd(x, y), len(x["contacts"] & y["contacts"]) / len(x["contacts"] | y["contacts"]))
    groups_by_threshold = {}
    for cutoff in (2, 3, 5):
        groups = [(cid,) for cid in sorted(poses)]
        while True:
            merges = []
            for i, j in combinations(range(len(groups)), 2):
                values = [pairs.get(tuple(sorted((a,b))), (float("inf"), 0))
                          for a in groups[i] for b in groups[j]]
                largest = max(d for d, _ in values)
                if largest <= cutoff and min(s for _, s in values) >= 0.5:
                    merges.append((largest, i, j))
            if not merges:
                break
            _, i, j = min(merges)
            groups[i] = tuple(sorted(groups[i] + groups.pop(j)))
        groups_by_threshold[str(cutoff)] = [
            dict(members=group, representative=min(group, key=lambda c: (-poses[c]["ipsae"], c)),
                 site=poses[group[0]]["site"]) for group in groups]
    return dict(job_id=os.environ["SLURM_JOB_ID"], cohort_size=len(poses), anchors=anchors,
                groups_by_rmsd_A=groups_by_threshold, contact_Jaccard_min=0.5, heavy_atom_contact_A=4.5,
                pairs=[dict(a=a,b=b,framework_rmsd_A=d,contact_jaccard=j) for (a,b),(d,j) in pairs.items()],
                poses={k: {**{f: v for f,v in p.items() if f not in ("anchor", "frame", "contacts")},
                           "canonical_contacts": sorted(p["contacts"])} for k,p in poses.items()})


def audit(row, mature):
    cid = row["design_id"]
    cp = ROOT / "predictions" / (cid + "_complex")
    mp = ROOT / "predictions" / (cid + "_monomer")
    s, m = read(cp / "summary.json"), read(mp / "summary.json")
    protocol = read(ROOT / "protocol.json")
    for record in (s, m):
        assert record["ok"] and record["pae_cif_mapping_validated"]
        assert record["protocol_id"] == protocol["protocol_id"]
        assert record["esm"] == "3.4.1.post1"
        assert not record["templates"] and not record["initial_guess"]
        assert all(protocol["predictor"][k] == v for k, v in record["kwargs"].items())
    assert [r["chain"] for r in s["msa_records"]] == ["A", "B"]
    assert s["msa_records"][0]["sha256"] == protocol["target"]["msa_sha256"]
    assert s["msa_records"][0]["rows"] > 1 and s["msa_records"][1]["rows"] == 0
    assert len(m["msa_records"]) == 1 and m["msa_records"][0]["rows"] == 0
    assert set(s["ipsae_directions"]) == {"A_to_B", "B_to_A"}
    assert s["ipsae_min"] == min(s["ipsae_directions"].values())
    atoms = pdb(cp / "model.pdb")
    t, tc = chain(atoms, "A", mature)
    b, bc = chain(atoms, "B", row["sequence"])
    renamed = atoms.copy()
    renamed.chain_id[atoms.chain_id == "A"] = "T"
    renamed.chain_id[atoms.chain_id == "B"] = "H"
    f = PDBFile()
    f.set_structure(renamed)
    f.write(cp / "mature_pose.pdb")
    dock = cp / "self_consistency_dockq.json"
    result = subprocess.run([str(DOCKQ), str(cp / "mature_pose.pdb"), row["reference_pdb"],
                             "--mapping", "TH:TH", "--n_cpu", "1", "--json", str(dock)],
                            capture_output=True, text=True, timeout=45)
    (cp / "dockq_stderr.txt").write_text(result.stderr)
    result.check_returncode()
    dq = read(dock)
    assert len(dq["best_result"]) == 1
    dockq = float(dq["GlobalDockQ"])
    mono, mc = chain(pdb(mp / "model.pdb"), "B", row["sequence"])
    match = frame_match(row)
    sg = mono[(mono.res_name == "CYS") & (mono.atom_name == "SG")]
    assert len(sg) == row["sequence"].count("C")
    d = distances(sg.coord, sg.coord)
    core_ids = [mc.res_id[len(row["framework_core"][0])-1], mc.res_id[match.start(3)-1]]
    i, j = [list(sg.res_id).index(r) for r in core_ids]

    def paired(indices):
        if not indices:
            return True
        first, *rest = indices
        return any(d[first, k] <= 3.0 and paired([r for r in rest if r != k]) for k in rest)

    sulfur_ok = bool(d[i,j] <= 3.0 and paired([k for k in range(len(sg)) if k not in (i,j)]))
    t, b = heavy(t), heavy(b)
    contact = distances(b.coord, t.coord).min(axis=1) <= 4.5
    interface = np.isin(bc.res_id, b.res_id[contact])
    cdr = np.zeros(len(bc), dtype=bool)
    for g in (1, 2, 3):
        cdr[match.start(g):match.end(g)] = True

    def anchor_distance(pos, atom):
        p = t[(t.res_id == tc.res_id[pos-30]) & (t.atom_name == atom)]
        assert len(p) == 1
        return float(distances(b.coord, p.coord).min())

    return dict(**row, v2_ipsae_min=s["ipsae_min"], v2_ipsae_directions=s["ipsae_directions"],
                delta_ipsae=s["ipsae_min"]-row["v1_ipsae_min"], v2_ipTM=s["ipTM"],
                v2_dockq=dockq, v2_monomer_plddt=m["chain_plddt"]["B"],
                monomer_disulfide_geometry_ok=sulfur_ok, monomer_core_SG_A=float(d[i,j]),
                cdr_contact_fraction=float(cdr[interface].mean()) if interface.any() else None,
                native_C_terminal_carbonyl_C_clearance_A=anchor_distance(211, "C"),
                glycan_anchor_ND2_clearance_A={p: anchor_distance(p, "ND2") for p in (54,145,151)},
                comparison_scope="target_and_VHH_S128_changed" if row["binder_changed"] else "target_input_changed_only",
                complex_cif_sha256=digest(cp / "model.cif"), monomer_cif_sha256=digest(mp / "model.cif"),
                prediction_summary_sha256=digest(cp / "summary.json"), design_reference_sha256=digest(row["reference_pdb"]),
                v2_quality_pass=(s["ipsae_min"] >= 0.6 and dockq >= 0.23
                                 and m["chain_plddt"]["B"] >= 70 and sulfur_ok), final_rank_eligible=False)


def main():
    assert os.environ.get("SLURM_JOB_ID"), "Slurm compute nodes only"
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=("cluster", "audit"))
    stage = p.parse_args().stage
    rows = read(ROOT / "inputs/cohort.json")
    canon = read(ROOT / "P19438.json")["sequence"]["value"]
    if stage == "cluster":
        out = ROOT / "reports/pose_clusters.json"
        if not out.exists():
            write(out, cluster(rows, canon))
    else:
        scored, failures = [], []
        for r in rows:
            try:
                scored.append(audit(r, canon[29:211]))
            except Exception as exc:
                failures.append(dict(id=r["design_id"], error=repr(exc)))
        write(ROOT / "reports/mature_core_scores.json",
              dict(job_id=os.environ["SLURM_JOB_ID"], scored=scored, failures=failures,
                   quality_pass_count=sum(r["v2_quality_pass"] for r in scored)))
        if failures:
            raise SystemExit(f"Incomplete audit: {len(failures)} failures; report retained")


if __name__ == "__main__":
    main()
