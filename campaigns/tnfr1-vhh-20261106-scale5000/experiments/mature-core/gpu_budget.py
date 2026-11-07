"""Fail closed unless accounting leaves room for the full two-hour validation."""
import json
import os
from pathlib import Path
import re
import subprocess

assert os.environ.get("SLURM_JOB_ID"), "Slurm compute nodes only"
main = Path(__file__).resolve().parents[2] / "tnfr1_vhh_20261106_scale5000"
conf = json.loads((main / "config/campaign.json").read_text())
rows = (main / "reports/submission_ids.tsv").read_text().splitlines()
ids = sorted({job for stage, job, wave in map(str.split, rows) if stage in {"generation", "scoring"}})
assert ids and all(i.isdecimal() for i in ids)
assert os.environ["SLURM_ARRAY_JOB_ID"] in ids
ledger = subprocess.check_output(["sacct", "-X", "-n", "-P", "-j", ",".join(ids),
                                  "--format=JobID%100,ElapsedRaw,AllocTRES"], text=True, timeout=30)
seen = set()
seconds = 0
for line in ledger.splitlines():
    job, elapsed, tres = [f.strip() for f in line.split("|")[:3]]
    if not re.fullmatch(r"\d+_\d+", job) or job in seen:
        continue
    gpus = re.search(r"(?:^|,)gres/gpu=(\d+)(?:,|$)", tres)
    if gpus:
        seen.add(job)
        seconds += int(elapsed) * int(gpus[1])
assert seconds > 0, "Campaign accounting missing"
recorded = json.loads((main / "reports/slurm_cost.json").read_text())["cumulative_allocated_gpu_hours"]
total = max(recorded, conf["baseline_allocated_gpu_hours"] + seconds / 3600)
print(f"GPU budget: cumulative {total:.3f} GPUh; this job reserves at most 2 additional GPUh", flush=True)
assert total + 2 <= conf["new_gpu_stage_stop_threshold"], "GPU budget guard"
