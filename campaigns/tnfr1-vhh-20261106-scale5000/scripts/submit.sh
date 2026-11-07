#!/bin/bash
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000
mkdir "$ROOT/reports/submission.lock"
IDS="$ROOT/reports/submission_ids.tsv"
: > "$IDS"
PREP=$(sbatch --parsable --hold "$ROOT/scripts/prepare.slurm")
printf 'prepare\t%s\t-1\n' "$PREP" >> "$IDS"
bash "$ROOT/scripts/submit_wave.sh" 0 "$PREP"
scontrol release "$PREP"
MON=$(sbatch --parsable "$ROOT/scripts/monitor.slurm")
printf 'monitor\t%s\t-1\n' "$MON" >> "$IDS"
cat "$IDS"
