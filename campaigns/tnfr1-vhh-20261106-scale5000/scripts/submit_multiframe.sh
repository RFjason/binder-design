#!/bin/bash
# Login-node Slurm management only. Caller must have inspected the CPU preflight.
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000
grep -q '"ok": true' "$ROOT/reports/multiframe_preflight.json"
mkdir "$ROOT/reports/multiframe_submission.lock"
IDS="$ROOT/reports/submission_ids.tsv"
# Combined primary+pilot limits remain generation6 / scoring2.
scontrol update JobId=64675207 ArrayTaskThrottle=4
scontrol update JobId=64675209 ArrayTaskThrottle=1
GEN=$(sbatch --parsable --hold --array=306-323%2 "$ROOT/scripts/generation.slurm")
printf 'generation\t%s\t-2\n' "$GEN" >> "$IDS"
COL=$(sbatch --parsable --hold --array=306-323%1 "$ROOT/scripts/collection.slurm")
printf 'collection\t%s\t-2\n' "$COL" >> "$IDS"
SCORE=$(sbatch --parsable --hold --array=306-323%1 "$ROOT/scripts/scoring.slurm")
printf 'scoring\t%s\t-2\n' "$SCORE" >> "$IDS"
AUDIT=$(sbatch --parsable --hold --array=306-323%2 "$ROOT/scripts/audit.slurm")
printf 'audit\t%s\t-2\n' "$AUDIT" >> "$IDS"
while IFS=$'\t' read -r task first; do
 DEP=''
 if [[ "$task" != "$first" ]]; then DEP="afterok:${AUDIT}_$first"; fi
 scontrol update "JobId=${GEN}_$task" Partition=a800 "Dependency=$DEP"
 scontrol update "JobId=${COL}_$task" "Dependency=afterany:${GEN}_$task"
 scontrol update "JobId=${SCORE}_$task" Partition=a800 "Dependency=afterok:${COL}_$task"
 scontrol update "JobId=${AUDIT}_$task" "Dependency=afterany:${SCORE}_$task,afterok:${COL}_$task"
done < "$ROOT/config/multiframe.tsv"
printf '{"generation":"%s","collection":"%s","scoring":"%s","audit":"%s","pilot_generation_concurrency":2,"pilot_scoring_concurrency":1}\n' "$GEN" "$COL" "$SCORE" "$AUDIT" > "$ROOT/reports/multiframe_submission.json"
for JOB in "$AUDIT" "$SCORE" "$COL" "$GEN"; do scontrol release "$JOB"; done
cat "$ROOT/reports/multiframe_submission.json"
