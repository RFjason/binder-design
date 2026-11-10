#!/bin/bash
# Login node: stage files and submit to Slurm; never execute analysis here.
set -euo pipefail
ROOT=${1:?experiment root required}
MAIN=$(dirname "$ROOT")/tnfr1_vhh_20261106_scale5000
cd "$ROOT"
test -s inputs/cohort.json
test -s inputs/manifest.json
test -s provenance.json
test -s "$MAIN/reports/submission_ids.tsv"
sha256sum -c --quiet stage.sha256
while IFS= read -r path; do
    test -s "$path" || { echo "Missing historical input: $path" >&2; exit 1; }
done < inputs/source_paths.txt
bash -n scripts/cluster.slurm
bash -n scripts/rescore.slurm
mkdir -p logs reports
mkdir reports/submission.lock
cpu=$(sbatch --parsable --hold --output="$ROOT/logs/cluster_%A_%a.out" scripts/cluster.slurm "$ROOT")
cpu=${cpu%%;*}
printf 'audit\t%s\t-4\n' "$cpu" | tee -a reports/submitted.tsv >> "$MAIN/reports/submission_ids.tsv"
gpu=$(sbatch --parsable --hold --output="$ROOT/logs/rescore_%A_%a.out" scripts/rescore.slurm "$ROOT")
gpu=${gpu%%;*}
printf 'scoring\t%s\t-4\n' "$gpu" | tee -a reports/submitted.tsv >> "$MAIN/reports/submission_ids.tsv"
# One slot moves from the current production generation array to validation.
scontrol update JobId=64675207 ArrayTaskThrottle=3
scontrol release "$cpu"
scontrol release "$gpu"
cat reports/submitted.tsv
squeue -h -r -j "$cpu,$gpu" -o '%i|%P|%T|%R'
