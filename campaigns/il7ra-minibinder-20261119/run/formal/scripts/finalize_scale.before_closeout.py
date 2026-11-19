"""Export the completed scoring pool and the current structural representatives."""
import csv
import json
import os
from pathlib import Path

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
root = Path(__file__).resolve().parents[1]
config = json.loads((root / "config/campaign.json").read_text())
candidates = json.loads((root / "reports/candidates.json").read_text())
audit = json.loads((root / "reports/pilot_audit.json").read_text())
shortlist = json.loads((root / "reports/combined_shortlist.json").read_text())
assert len(candidates) == audit["scored"] == config["new_quota"]
assert len({row["sequence"] for row in candidates}) == len(candidates)
assert config["baseline_count"] + audit["scored"] == config["target_count"]
fields = ["id", "sequence", "ipSAE_min", "monomer_pLDDT", "crystal_glycan_min_A",
          "fold_cluster", "interface_cluster", "ipSAE_min_M1_G236"]
for row in shortlist:
    path = root / "predictions" / (row["id"] + "_historical") / "summary.json"
    row["ipSAE_min_M1_G236"] = json.loads(path.read_text())["ipsae_min"] if path.is_file() else None
with (root / "reports/shortlist.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(shortlist)
(root / "reports/shortlist.fasta").write_text("".join(
    f">{row['id']}\n{row['sequence']}\n" for row in shortlist))
with (root / "reports/scored_candidates.jsonl").open("w") as handle:
    for row in audit["rows"]:
        handle.write(json.dumps(row) + "\n")
diversity = json.loads((root / "reports/diversity.json").read_text())
summary = {"total_scored": config["target_count"], "baseline": config["baseline_count"],
           "new_scored": audit["scored"], "passing_new": sum(r["pilot_gate"] for r in audit["rows"]),
           "fold_clusters_passing_combined": diversity["fold_clusters"],
           "interface_clusters_passing_combined": diversity["interface_clusters"],
           "representatives": len(shortlist), "historical_rescored": sum(
               r["ipSAE_min_M1_G236"] is not None for r in shortlist),
           "remaining_checks": ["Full glycan accessibility", "Complete target-specific patent and sequence novelty"]}
(root / "reports/completed.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
