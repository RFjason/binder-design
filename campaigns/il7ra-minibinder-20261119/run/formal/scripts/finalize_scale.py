"""Export the completed scoring pool and the current structural representatives."""
import csv
import fcntl
import json
import os
from pathlib import Path

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
root = Path(__file__).resolve().parents[1]
config = json.loads((root / "config/campaign.json").read_text())
candidates = json.loads((root / "reports/candidates.json").read_text())
audit = json.loads((root / "reports/pilot_audit.json").read_text())
shortlist = json.loads((root / "reports/combined_shortlist.json").read_text())
closeout_path = root / "state/closeout.json"
closeout = json.loads(closeout_path.read_text()) if closeout_path.is_file() else None
expected = closeout["new_scored_target"] if closeout else config["new_quota"]
assert len(candidates) == audit["scored"] == expected
assert len({row["sequence"] for row in candidates}) == len(candidates)
assert {row["id"] for row in candidates} == {row["id"] for row in audit["rows"]}
baseline = json.loads(Path(config["baseline_audit"]).read_text())["rows"]
assert len(baseline) == config["baseline_count"]
scored = baseline + audit["rows"]
assert len({row["sequence"] for row in scored}) == len(scored) == config["baseline_count"] + expected
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
    for row in scored:
        handle.write(json.dumps(row) + "\n")
score_fields = ["id", "sequence", "ipSAE_min", "monomer_pLDDT", "monomer_rmsd_A",
                "crystal_glycan_min_A", "sequon_ND2_min_A", "unresolved_end_min_A",
                "interchain_clash_pairs", "paired_cysteines", "pilot_gate", "shortlisted"]
shortlist_ids = {row["id"] for row in shortlist}
with (root / "reports/scored_candidates.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, score_fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows({**row, "shortlisted": row["id"] in shortlist_ids}
                     for row in sorted(scored, key=lambda row: -row["ipSAE_min"]))
diversity = json.loads((root / "reports/diversity.json").read_text())
summary = {"total_scored": len(scored), "original_target": config["target_count"],
           "generation_stopped_early": closeout is not None, "baseline": config["baseline_count"],
           "new_scored": audit["scored"], "passing_new": sum(r["pilot_gate"] for r in audit["rows"]),
           "passing_total": sum(r["pilot_gate"] for r in scored), "audit_failures": len(audit["failures"]),
           "fold_clusters_passing_combined": diversity["fold_clusters"],
           "interface_clusters_passing_combined": diversity["interface_clusters"],
           "representatives": len(shortlist), "historical_rescored": sum(
               r["ipSAE_min_M1_G236"] is not None for r in shortlist),
           "remaining_checks": ["Full glycan accessibility", "Complete target-specific patent and sequence novelty"]}
(root / "reports/completed.json").write_text(json.dumps(summary, indent=2) + "\n")
state_path = root / "state/controller.json"
with (root / "state/controller.lock").open("a") as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    state = json.loads(state_path.read_text())
    state["status"] = "CLOSEOUT_COMPLETE" if closeout else "COMPLETE"
    temporary = state_path.with_suffix(".json.finalize.tmp")
    temporary.write_text(json.dumps(state, indent=2) + "\n")
    temporary.replace(state_path)
print(json.dumps(summary), flush=True)
