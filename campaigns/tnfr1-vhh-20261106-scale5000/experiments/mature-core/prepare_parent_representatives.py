"""Stage unchanged scale5000 parent representatives under the frozen v2 protocol."""
import argparse
import json
from pathlib import Path
import re

from prepare import AA, FRAME, HERE, MAIN, MAIN_REMOTE, REMOTE, dump, jsonl, sha

ROOT = REMOTE / "tnfr1_parent_rescore_20261110"
SNAPSHOT = MAIN / "results/completed_20261110/scored_candidates.jsonl"
PRIOR_SCORES = HERE / "results/mature_core_scores.json"
PROTOCOL = HERE / "submission/protocol.json"
PRIOR_HASHES = HERE / "submission/provenance.json"
PRIORITY = (
    "bc2_gpu_crd2_08:trajectory_1",
    "rfd4_crd2_compact_00:tnfr1_crd2_compact_rep1_generated",
)


def build(prior, stage):
    stage.mkdir()
    (stage / "scripts").mkdir()
    (stage / "inputs").mkdir()
    pinned = json.loads(PRIOR_HASHES.read_bytes())["staged_sha256"]
    protocol = json.loads(PROTOCOL.read_bytes())
    msa_rel = protocol["target"]["msa"]
    for rel in ("P19438.json", msa_rel, "scripts/predict_esmfold2.py",
                "scripts/analyze.py", "scripts/gpu_budget.py"):
        data = (prior / rel).read_bytes()
        assert sha(data) == pinned[rel], rel
        (stage / rel).write_bytes(data)
    canonical = json.loads((stage / "P19438.json").read_bytes())["sequence"]["value"]
    target = canonical[29:211]
    assert sha(target.encode()) == protocol["target"]["mature_sequence_sha256"]
    assert sha((stage / msa_rel).read_bytes()) == protocol["target"]["msa_sha256"]
    protocol["remote_root"] = str(ROOT)
    dump(stage / "protocol.json", protocol)
    wrapper = (prior / "scripts/rescore.slurm").read_bytes()
    assert sha(wrapper) == pinned["scripts/rescore.slurm"]
    wrapper = wrapper.decode().replace("tnfr1_mature21", "tnfr1_parent18")
    wrapper = wrapper.replace("02:00:00", "00:45:00").replace("98m", "32m")
    wrapper = "\n".join(line for line in wrapper.split("\n")
                         if not line.startswith(("# If the CPU partition", "flock ")))
    (stage / "scripts/rescore.slurm").write_bytes(wrapper.encode())

    prior_scores = json.loads(PRIOR_SCORES.read_bytes())["scored"]
    prior_by_id = {r["source_design_id"]: r for r in prior_scores if not r["binder_changed"]}
    passing = [r for r in jsonl(SNAPSHOT) if r["computational_quality_pass"]]
    assert len(passing) == 97
    ordered = sorted(passing, key=lambda r: (-r["ipsae_min"], r["design_id"]))
    best = {}
    for r in ordered:
        best.setdefault(r["parent_id"], r)
    assert len(best) == 20 and all(p in best for p in PRIORITY)
    parents = list(PRIORITY) + [p for p in best if p not in PRIORITY]
    rows, cases, reused = [], [], []
    for parent in parents:
        r = best[parent]
        cid, seq = r["design_id"], r["sequence"]
        frame = (*FRAME[:3], FRAME[3] + "S")
        assert r["framework_id"] == "h-NbBCII10" and list(frame) == r["framework_core"]
        assert re.fullmatch(f"([{AA}]{{5,30}})".join(map(re.escape, frame)), seq)
        assert re.fullmatch(r"TNFR1_[A-Za-z0-9_]+", cid) and 100 <= len(seq) <= 150
        key = sha(seq.encode())
        if cid in prior_by_id:
            old = prior_by_id[cid]
            assert old["sequence_sha256"] == key and old["parent_id"] == parent
            reused.append(dict(source_design_id=cid, parent_id=parent, sequence_sha256=key,
                               v2_ipsae_min=old["v2_ipsae_min"], v2_quality_pass=old["v2_quality_pass"]))
            continue
        variant = cid + "_mature_core"
        reference = MAIN_REMOTE / r["reference_pdb"]
        assert reference.is_relative_to(MAIN_REMOTE / "candidates")
        rows.append(dict(design_id=variant, source_design_id=cid, sequence=seq, length=len(seq),
                         source_sequence=seq, binder_changed=False, sequence_sha256=key,
                         source_sequence_sha256=key, framework_core=frame, reference_pdb=str(reference),
                         stratum=r["stratum"], method=r["structure_method"], seq_method=r["seq_method"],
                         parent_id=parent, v1_ipsae_min=r["ipsae_min"], v1_pose_dockq=r["pose_dockq"],
                         v1_monomer_plddt=r["monomer_plddt"], source_generation_sha256=r["source_sha256"]))
        binder = dict(id="B", sequence=seq)
        chains = [dict(id="A", sequence=target, msa=msa_rel), binder]
        cases += [dict(id=variant + "_complex", candidate_id=variant, chains=chains),
                  dict(id=variant + "_monomer", candidate_id=variant, chains=[binder])]
    assert len(rows) == len({r["sequence_sha256"] for r in rows}) == 18
    assert len(reused) == 2 and all(r["v2_quality_pass"] for r in reused)
    assert tuple(r["parent_id"] for r in rows[:2]) == PRIORITY
    dump(stage / "inputs/cohort.json", rows)
    dump(stage / "inputs/manifest.json", cases)
    dump(stage / "inputs/reused.json", reused)
    (stage / "inputs/reference_paths.txt").write_bytes(
        ("\n".join(r["reference_pdb"] for r in rows) + "\n").encode())
    files = {p.relative_to(stage).as_posix(): sha(p.read_bytes())
             for p in sorted(stage.rglob("*")) if p.is_file()}
    dump(stage / "provenance.json", {
        "selection": "highest v1 ipSAE_min per passing scale5000 parent; exact prior v2 matches reused",
        "source_sha256": {p.relative_to(MAIN).as_posix(): sha(p.read_bytes())
                          for p in (SNAPSHOT, PRIOR_SCORES, PROTOCOL, PRIOR_HASHES)},
        "staged_sha256": files,
    })
    (stage / "stage.sha256").write_bytes("".join(f"{h}  {p}\n" for p, h in files.items()).encode())
    print(f"Validated {len(rows)} new representatives, {len(reused)} reused and {len(cases)} prediction cases")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("prior_stage", type=Path)
    parser.add_argument("stage", type=Path)
    args = parser.parse_args()
    build(args.prior_stage, args.stage)
