#!/bin/bash
#SBATCH --partition=64c512g
#SBATCH --exclude=node019
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000/logs/grouped_cpu_%A_%a.out
set -euo pipefail
cd "${1:?campaign root required}"
stage=${2:?stage required}
count=${3:?group size required}
[[ "$stage" == collection || "$stage" == audit ]]
[[ "$count" =~ ^[1-9][0-9]*$ ]]
start=${SLURM_ARRAY_TASK_ID:?Slurm array task required}
failed=0
while IFS=$'\t' read -r task method batch genpart scorepart canary; do
  [[ "$task" =~ ^[0-9]+$ ]] || continue
  (( task >= start && task < start+count && task <= 305 )) || continue
  if [[ "$stage" == audit ]] && ! grep -Eq '"state":[[:space:]]*"COMPLETED"' "reports/jobs/${batch}_collection.json"; then
    echo "SKIP_AUDIT $task $batch collection incomplete"
    continue
  fi
  if ! S5K_TASK="$task" bash "scripts/$stage.slurm"; then
    echo "GROUP_MEMBER_FAILED $stage $task $batch"
    failed=1
  fi
done < config/batches.tsv
# Individual records drive scoring and auditing. One collection failure
# must not block the other members of this group.
if [[ "$stage" == audit ]]; then exit "$failed"; fi
