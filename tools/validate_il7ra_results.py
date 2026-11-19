"""Validate IL7RA published evidence with the Python standard library only."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import zipfile

REPO = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--campaign", type=Path, default=REPO / "campaigns/il7ra-minibinder-20261119")
parser.add_argument("--archives", type=Path)
args = parser.parse_args()
root = args.campaign.resolve()


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fasta(path):
    records = {}
    key = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            assert key not in records
            records[key] = ""
        elif line.strip():
            assert key is not None
            records[key] += line.strip()
    return records


AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
              "ARNDCQEGHILKMFPSTWYV"))


def pdb_sequences(path):
    chains = {}
    seen = set()
    for line in path.read_text().splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
            key = line[21:27]
            if key not in seen:
                seen.add(key)
                chain = line[21]
                chains[chain] = chains.get(chain, "") + AA[line[17:20]]
    return chains


manifest = load(root / "provenance/data_manifest.json")
for name, record in manifest["files"].items():
    path = (root / name).resolve()
    assert path.is_relative_to(root) and path.is_file(), name
    assert path.stat().st_size == record["bytes"] and sha(path) == record["sha256"], name
summary = load(root / "results/completed.json")
rows = [json.loads(line) for line in (root / "results/scored_candidates.jsonl").read_text().splitlines()]
assert len(rows) == len({r["id"] for r in rows}) == len({r["sequence"] for r in rows}) == 2855
rowmap = {r["id"]: r for r in rows}
passing = {r["id"]: r for r in rows if r["pilot_gate"]}
assert len(passing) == summary["passing_total"] == 445
for row in rows:
    assert 55 <= len(row["sequence"]) <= 95 and set(row["sequence"]) <= set(AA.values())
    assert hashlib.sha256(row["sequence"].encode()).hexdigest()[:16] == row["id"].rsplit("_", 1)[-1]
    assert abs(min(row["ipSAE_directions"].values()) - row["ipSAE_min"]) < 1e-6
    expected = (row["ipSAE_min"] >= .6 and row["monomer_pLDDT"] >= 75 and row["monomer_rmsd_A"] <= 4
                and row["crystal_glycan_min_A"] >= 3.2 and row["sequon_ND2_min_A"] >= 4.5
                and row["unresolved_end_min_A"] >= 3.2 and row["interchain_clash_pairs"] == 0
                and row["paired_cysteines"] and bool(row["passing_sources"])
                and not row["initial_public_reference_hits"])
    assert row["pilot_gate"] == expected, row["id"]
    for source in row["passing_sources"]:
        assert source["pose_ca_rmsd_A"] <= 5 and source["design_contact_recall"] >= .4
        assert len(source["hotspots_contacted"]) >= source.get("min_hotspots", 2)
with (root / "results/scored_candidates.csv").open(newline="") as f:
    csvrows = list(csv.DictReader(f))
assert len(csvrows) == 2855 and {r["id"] for r in csvrows} == set(rowmap)
with (root / "passing/candidates.csv").open(newline="") as f:
    exported = list(csv.DictReader(f))
assert len(exported) == 445 and {r["id"] for r in exported} == set(passing)
assert fasta(root / "passing/candidates.fasta") == {cid: r["sequence"] for cid, r in passing.items()}
diversity = load(root / "results/diversity.json")
assert {r["id"] for r in diversity["members"]} == set(passing)
assert diversity["fold_clusters"] == 209 and diversity["interface_clusters"] == 18
chosen, parents, pairs = [], set(), set()
for member in sorted(diversity["members"], key=lambda x: -x["ipSAE_min"]):
    source = next((s for s in member["sources"] if s["parent"] not in parents), None)
    pair = member["fold_cluster"], member["interface_cluster"]
    if source and pair not in pairs and len(chosen) < 20:
        chosen.append(member["id"])
        parents.add(source["parent"])
        pairs.add(pair)
selected = load(root / "top20/selection.json")["selected"]
assert [r["id"] for r in selected] == chosen
assert len(chosen) == len(parents) == len(pairs) == 20
assert fasta(root / "top20/candidates.fasta") == {r["id"]: r["sequence"] for r in selected}
protocol = load(root / "run/formal/scoring_protocol.json")
target = fasta(root / "run/formal/inputs/IL7RA_E21_D239.fasta")
assert len(target) == 1
target = next(iter(target.values()))
assert len(target) == 219 and hashlib.sha256(target.encode()).hexdigest() == protocol["target"]["sequence_sha256"]
assert sha(root / "run/formal/inputs/IL7RA_E21_D239.a3m") == protocol["target"]["msa"]["sha256"]
assert protocol["predictor"]["seed"] == 0 and protocol["predictor"]["num_diffusion_samples"] == 1
for cid, row in passing.items():
    for kind in ("complex", "monomer"):
        folder = root / "passing/structures" / cid / kind
        pred = load(folder / "summary.json")
        assert pred["ok"] and pred["pae_cif_mapping_validated"]
        assert pred["kwargs"]["seed"] == 0 and not pred["templates"] and not pred["initial_guess"]
        chains = pdb_sequences(folder / "model.pdb")
        assert chains == ({"A": target, "B": row["sequence"]} if kind == "complex" else {"B": row["sequence"]}), cid
        if kind == "complex":
            assert abs(pred["ipsae_min"] - row["ipSAE_min"]) < 1e-6
        else:
            assert abs(pred["chain_plddt"]["B"] - row["monomer_pLDDT"]) < 1e-5
for row in selected:
    result = load(root / "top20/historical" / row["id"] / "summary.json")
    assert result["ok"] and abs(result["ipsae_min"] - row["ipSAE_min_M1_G236"]) < 1e-6
    assert row["ipSAE_min_M1_G236"] >= .6
provenance = load(root / "provenance/passing_sources.json")
assert len(provenance) == 445 and {r["id"] for r in provenance} == set(passing)
if args.archives:
    for part in load(root / "provenance/archive_manifest.json")["parts"]:
        archive = args.archives / part["filename"]
        assert archive.stat().st_size == part["bytes"] and sha(archive) == part["sha256"]
        with zipfile.ZipFile(archive) as handle:
            assert handle.testzip() is None
            names = set(handle.namelist())
            assert len(names) == part["entries"]
            for row in provenance:
                for source in row["source_files"]:
                    assert hashlib.sha256(handle.read(source["archive_path"])).hexdigest() == source["sha256"]
                for kind in ("complex", "monomer"):
                    required = ["model.pdb", "model.cif", "summary.json", "pae.npz", "input_and_token_mapping.json"]
                    if kind == "complex":
                        required += ["af3like_pae.json", "model_10_10.txt", "model_10_10_byres.txt"]
                    for filename in required:
                        assert f"evidence/{row['id']}/{kind}/{filename}" in names
print(f"PASS IL7RA: 2855 unique scored sequences; 445 passing; 20 deterministic representatives; {len(manifest['files'])} files verified.")
