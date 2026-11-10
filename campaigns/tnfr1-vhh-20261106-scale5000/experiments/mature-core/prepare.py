"""Build the frozen 21-candidate validation package from the public snapshots."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re

HERE = Path(__file__).resolve().parent
MAIN = HERE.parents[1]
OLD = MAIN.parent / "tnfr1-vhh-20261103"
REMOTE = PurePosixPath("/dssg/home/acct-clshzh/jasoncao1280/benchmarks")
OLD_REMOTE = REMOTE / "tnfr1_vhh_20261103"
MAIN_REMOTE = REMOTE / "tnfr1_vhh_20261106_scale5000"
ROOT = REMOTE / "tnfr1_mature_rescore_20261107"
FRAME = ("QVQLVESGGGLVQPGGSLRLSC", "WFRQAPGQGLEAVA",
         "YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC", "WGQGTLVTVS")
AA = "ACDEFGHIKLMNPQRSTVWY"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def jsonl(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def dump(path, data):
    path.write_bytes((json.dumps(data, indent=2, ensure_ascii=True) + "\n").encode())


def build(msa_input, stage):
    inputs = stage / "inputs"
    scripts = stage / "scripts"
    inputs.mkdir(parents=True, exist_ok=True)
    scripts.mkdir(exist_ok=True)
    canonical_json = (OLD / "inputs/uniprot_P19438.json").read_bytes()
    canonical = json.loads(canonical_json)["sequence"]["value"]
    mature = canonical[29:211]
    assert len(mature) == 182 and mature[0] == "L" and mature[-1] == "T"
    msa = msa_input.read_bytes()
    seqs = [re.sub("[a-z.]", "", "".join(row.splitlines()[1:])) for row in msa.decode().split(">")[1:]]
    assert len(seqs) > 1 and seqs[0] == mature and {len(seq) for seq in seqs} == {182}
    (stage / "P19438.json").write_bytes(canonical_json)
    msa_rel = "inputs/TNFR1_mature_L30_T211.a3m"
    (stage / msa_rel).write_bytes(msa)
    protocol = json.loads((OLD / "config/scoring_protocol_v1.json").read_text())
    protocol.update(protocol_id="TNFR1_VHH_ESMFold2_L30_T211_core_v2_20261107", remote_root=str(ROOT))
    protocol["target"] = dict(uniprot="P19438-1", mature_range=[30, 211],
                              mature_sequence_sha256=sha(mature.encode()), msa=msa_rel,
                              msa_sha256=sha(msa), msa_rows=len(seqs), reference="9S3I chain A",
                              his_tag="Declared at target C terminus; sequence/linker unknown; omitted from this experiment")
    protocol["binder"]["framework"] = "h-NbBCII10, complete FR4 including IMGT S128"
    protocol["score"]["pose_check"] = "DockQ vs original mature target design pose, chain mapping TH:TH"
    dump(stage / "protocol.json", protocol)
    sources = [OLD / "results/all_candidates.jsonl",
               MAIN / "results/progress_20261107/scored_candidates.jsonl"]
    rows, cases = [], []
    for source, remote, expected_count in zip(sources, (OLD_REMOTE, MAIN_REMOTE), (19, 2)):
        selected = [r for r in jsonl(source) if r.get("computational_quality_pass") is True]
        assert len(selected) == expected_count
        for r in sorted(selected, key=lambda x: x["design_id"]):
            cid = r["design_id"]
            assert re.fullmatch(r"TNFR1_[A-Za-z0-9_]+", cid)
            old = r["sequence"]
            assert old.endswith("WGQGTLVTVS" if expected_count == 19 else "WGQGTLVTVSS")
            seq = old + "S" if expected_count == 19 else old
            assert 100 <= len(seq) <= 150 and set(seq) <= set(AA)
            frame = (*FRAME[:3], FRAME[3] + "S")
            assert re.fullmatch(f"([{AA}]{{5,30}})".join(map(re.escape, frame)), seq)
            variant = cid + "_mature_core"
            folder = remote / "predictions" / (cid + "_complex")
            row = dict(design_id=variant, source_design_id=cid, sequence=seq, length=len(seq),
                       source_sequence=old, binder_changed=(old != seq), sequence_sha256=sha(seq.encode()),
                       source_sequence_sha256=sha(old.encode()), framework_core=frame,
                       source_complex=str(folder / "model.cif"), source_model=str(folder / "model.pdb"),
                       reference_pdb=str(remote / r["reference_pdb"]),
                       stratum=r["stratum"], method=r["structure_method"], parent_id=r["parent_id"],
                       v1_ipsae_min=r["ipsae_min"], v1_pose_dockq=r["pose_dockq"],
                       v1_monomer_plddt=r["monomer_plddt"], source_generation_sha256=r["source_sha256"])
            rows.append(row)
            binder = dict(id="B", sequence=seq)
            target = dict(id="A", sequence=mature, msa=msa_rel)
            cases.extend(dict(id=variant + "_" + kind, candidate_id=variant, chains=chains)
                         for kind, chains in (("complex", [target, binder]), ("monomer", [binder])))
    assert len({r["sequence"] for r in rows}) == len(rows) == 21
    dump(inputs / "cohort.json", rows)
    dump(inputs / "manifest.json", cases)
    paths = sorted({r[key] for r in rows for key in ("source_complex", "source_model", "reference_pdb")})
    assert all(re.fullmatch(r"/[A-Za-z0-9_./-]+", p) for p in paths)
    (inputs / "source_paths.txt").write_bytes(("\n".join(paths) + "\n").encode())
    predictor_path = MAIN / "scripts/predict_esmfold2.py"
    predictor = predictor_path.read_text()
    for before, after in (
        (f"Path('{MAIN_REMOTE}')", "Path(__file__).resolve().parents[1]"),
        ("'scoring_protocol_v1.json'", "'protocol.json'"),
        ("canonical[:211]", "canonical[29:211]"),
        ("'competition_sequence_sha256'", "'mature_sequence_sha256'"),
    ):
        assert predictor.count(before) == 1
        predictor = predictor.replace(before, after)
    (scripts / "predict_esmfold2.py").write_bytes(predictor.encode())
    for name in ("analyze.py", "gpu_budget.py", "cluster.slurm", "rescore.slurm", "submit.sh"):
        (scripts / name).write_bytes((HERE / name).read_text(encoding="utf-8").encode())
    dump(stage / "provenance.json", {
        "source_sha256": {p.relative_to(MAIN.parent).as_posix(): sha(p.read_bytes())
                          for p in [*sources, predictor_path, OLD / "config/scoring_protocol_v1.json"]},
        "staged_sha256": {p.relative_to(stage).as_posix(): sha(p.read_bytes())
                          for p in sorted(stage.rglob("*"))
                          if p.is_file() and p.suffix != ".pyc"
                          and p.name not in {"provenance.json", "stage.sha256"}},
    })
    manifest = json.loads((stage / "provenance.json").read_bytes())["staged_sha256"]
    (stage / "stage.sha256").write_bytes("".join(f"{h}  {p}\n" for p,h in manifest.items()).encode())
    print(f"Validated {len(rows)} distinct VHHs / {len(cases)} prediction cases / {len(seqs)} aligned target MSA rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("target_msa", type=Path)
    parser.add_argument("stage", type=Path)
    args = parser.parse_args()
    build(args.target_msa, args.stage)
