"""Preflight input hashes and validate the MSA with the actual ESM parser."""
import ast
import hashlib
import json
import os
from pathlib import Path

from esm.utils.msa import MSA

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
root = Path(__file__).resolve().parents[1]
protocol = json.loads((root / "scoring_protocol.json").read_text())
target = protocol["target"]
sequence = (root / target["fasta"]).read_text().splitlines()[1]
msa_path = root / target["msa"]["file"]
assert hashlib.sha256(sequence.encode()).hexdigest() == target["sequence_sha256"]
assert hashlib.sha256(msa_path.read_bytes()).hexdigest() == target["msa"]["sha256"]
assert len(sequence) == 219 and (sequence[0], sequence[45], sequence[-1]) == ("E", "I", "D")
msa = MSA.from_a3m(msa_path, remove_insertions=True)
assert msa.sequences[0] == sequence and msa.depth > 1
for path in sorted((root / "scripts").glob("*.py")):
    ast.parse(path.read_text(), filename=str(path))
controls = [case for i in range(3) for case in json.loads((root / "inputs" / f"controls_{i}.json").read_text())]
assert len({case["id"] for case in controls}) == 12
assert all(case["control_only"] and case["chains"][0]["sequence"] == sequence
           and "msa" not in case["chains"][1] for case in controls)
report = {"ok": True, "job": os.environ["SLURM_JOB_ID"], "length": len(sequence),
          "target_msa_rows_parsed": msa.depth, "controls": len(controls)}
(root / "reports/input_preflight.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report), flush=True)
