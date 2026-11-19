"""Export the prescreen and separate old-construct score for each selected parent."""
import csv
import json
import os
from pathlib import Path

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
ROOT = Path(__file__).resolve().parents[1]
reports = ROOT / "reports"
audit = json.loads((reports / "pilot_audit.json").read_text())
fields = ["id", "length", "method", "patch", "parent", "ipSAE_min",
          "ipSAE_min_M1_G236", "monomer_pLDDT", "crystal_glycan_min_A",
          "hotspots_contacted", "sequence"]
out = []
for row in audit["independent_shortlist"]:
    source = row["selected_source"]
    path = ROOT / "predictions" / (row["id"] + "_historical") / "summary.json"
    old = json.loads(path.read_text())["ipsae_min"] if path.is_file() else None
    out.append({"id": row["id"], "length": len(row["sequence"]),
                "method": source["method"], "patch": source["patch"], "parent": source["parent"],
                "ipSAE_min": row["ipSAE_min"], "ipSAE_min_M1_G236": old,
                "monomer_pLDDT": row["monomer_pLDDT"], "crystal_glycan_min_A": row["crystal_glycan_min_A"],
                "hotspots_contacted": ",".join(map(str, source["hotspots_contacted"])),
                "sequence": row["sequence"]})
with (reports / "pilot_shortlist.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(out)
(reports / "pilot_shortlist.fasta").write_text("".join(f">{row['id']}\n{row['sequence']}\n" for row in out))
summary = {"scored": audit["scored"], "audit_failures": len(audit["failures"]),
           "pilot_gate": sum(row["pilot_gate"] for row in audit["rows"]), "independent_shortlist": len(out),
           "historical_rescored": sum(row["ipSAE_min_M1_G236"] is not None for row in out),
           "controls_scored": sum(row["ipSAE_min"] is not None for row in audit["controls"]),
           "scope": "Computational pilot. Full glycans and patent novelty pending."}
(reports / "pilot_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
