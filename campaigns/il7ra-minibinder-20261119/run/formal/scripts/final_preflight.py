"""Verify the isolated TM-align wheel and all campaign entry points on a CPU node."""
import ast
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
from biotite.structure.io.pdb import PDBFile

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
root = Path(__file__).resolve().parents[1]
wheel = next((root / "sources").glob("tmtools-0.3.0-cp312-*.whl"))
assert hashlib.sha256(wheel.read_bytes()).hexdigest() == "74bdea863afa4c6c29b95635fa63c6b3f8eac32c996b95e22be1708c5e4e1205"
vendor = root / ".vendor"
vendor.mkdir(exist_ok=True)
with zipfile.ZipFile(wheel) as archive:
    assert all(not Path(name).is_absolute() and ".." not in Path(name).parts for name in archive.namelist())
    archive.extractall(vendor)
sys.path.insert(0, str(vendor))
from tmtools import tm_align

atoms = PDBFile.read(root / "inputs/receptor.pdb").get_structure(model=1)
ca = atoms[(atoms.chain_id == "A") & (atoms.atom_name == "CA")]
sequence = (root / "inputs/IL7RA_E21_D239.fasta").read_text().splitlines()[1]
coords = np.asarray(ca.coord, dtype=np.float64)
result = tm_align(coords, coords, sequence, sequence)
assert min(result.tm_norm_chain1, result.tm_norm_chain2) > 0.999
scripts = sorted((root / "scripts").glob("*.py"))
for path in scripts:
    ast.parse(path.read_text(), filename=str(path))
batches = json.loads((root / "inputs/design_batches.json").read_text())
counts = Counter()
for row in batches:
    path = root / "inputs" / (row["id"] + ".json")
    spec = json.loads(path.read_text())
    counts[row["method"], row["patch"].split("_")[0]] += row["replicates"]
    if row["method"] == "rfd4":
        assert f"A1-219" in spec[row["id"]]["contig"]
    else:
        assert spec["campaign_seed"] == row["seed"] and spec["max_trajectories"] == 1
assert [counts["rfd4", group] for group in ("known", "adjacent", "other")] == [2000, 1200, 800]
assert [counts["bc2", group] for group in ("known", "adjacent", "other")] == [150, 90, 60]
(root / "reports/pilot_audit.json").write_text(json.dumps({"rows": [], "scored": 0}))
env = {**os.environ, "SCALE_WAVE": "-1"}
for script in ("cluster_scale.py", "select_scale.py"):
    subprocess.run([sys.executable, str(root / "scripts" / script)], env=env, check=True)
report = {"ok": True, "scripts_parsed": len(scripts), "tm_align_identity": result.tm_norm_chain1,
          "batch_counts": {"/".join(k): v for k, v in counts.items()}, "job": os.environ["SLURM_JOB_ID"]}
(root / "reports/final_preflight.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report), flush=True)
