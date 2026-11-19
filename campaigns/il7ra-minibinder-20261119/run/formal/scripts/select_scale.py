"""Choose up to twenty representatives from combined fold and interface clusters."""
import json
import os
from pathlib import Path

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
root = Path(__file__).resolve().parents[1]
config = json.loads((root / "config/campaign.json").read_text())
new_audit = json.loads((root / "reports/pilot_audit.json").read_text())
baseline = json.loads(Path(config["baseline_audit"]).read_text())
diversity = json.loads((root / "reports/diversity.json").read_text())
rows = {r["id"]: r for r in baseline["rows"] + new_audit["rows"]}
ranked = sorted(diversity["members"], key=lambda member: -member["ipSAE_min"])
selected, parents, joints = [], set(), set()
for member in ranked:
    source = next((s for s in member["sources"] if s["parent"] not in parents), None)
    joint = member["fold_cluster"], member["interface_cluster"]
    if source and joint not in joints and len(selected) < 20:
        selected.append({**rows[member["id"]], **member, "selected_source": source})
        parents.add(source["parent"])
        joints.add(joint)
(root / "reports/combined_shortlist.json").write_text(json.dumps(selected, indent=2) + "\n")
old_sequence = "".join(line.strip() for line in (root / "inputs/IL7RA_competition_M1_G236.fasta").read_text().splitlines()
                       if not line.startswith(">"))
assert len(old_sequence) == 236
for shard in range(4):
    cases = [{"id": row["id"] + "_historical", "historical_recheck": True, "chains": [
        {"id": "A", "sequence": old_sequence, "msa": "inputs/IL7RA_competition_M1_G236.a3m"},
        {"id": "B", "sequence": row["sequence"]}]} for row in selected[shard::4]]
    (root / "inputs" / f"historical_score_{shard}.json").write_text(json.dumps(cases, indent=2) + "\n")
print(f"COMBINED_SHORTLIST {len(selected)} fold/interface representatives", flush=True)
