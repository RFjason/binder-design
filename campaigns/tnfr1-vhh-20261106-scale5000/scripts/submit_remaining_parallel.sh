#!/bin/bash
# Slurm administration only. Recover the held generation array from the
# first attempt, grouping CPU work to fit the 500-job association limit.
set -euo pipefail
ROOT=${1:?campaign root required}
GEN=${2:?held generation array required}
OLD_COL=${3:?held collection array required}
cd "$ROOT"
IDS=reports/submission_ids.tsv
GROUP_SIZE=8
[[ "$GEN" =~ ^[0-9]+$ && "$OLD_COL" =~ ^[0-9]+$ && "$GEN" != "$OLD_COL" ]]
for wave in 1 2 3; do test -d "reports/wave_${wave}.lock"; done
if awk -F '\t' '$3 > 0 {found=1} END {exit !found}' "$IDS"; then
  echo 'A later wave is already registered; refusing to submit duplicates' >&2
  exit 1
fi
grep -Fxq $'generation\t'"$GEN"$'\t3' reports/parallel_submission_attempt.tsv
grep -Fxq $'collection\t'"$OLD_COL"$'\t3' reports/parallel_submission_attempt.tsv
for entry in "$GEN:s5k_generation" "$OLD_COL:s5k_collection"; do
  IFS=: read -r id name <<< "$entry"
  [[ "$(squeue -h -j "$id" -o '%F|%j|%T|%R')" == "$id|$name|PENDING|(JobHeldUser)" ]]
done
gen_info=$(scontrol show job -o "$GEN")
grep -q 'ArrayTaskThrottle=4 ' <<< "$gen_info"
for route in rfd4 bc2_v2; do
  grep -Eq '"ok"[[:space:]]*:[[:space:]]*true' "reports/route_ready_$route.json"
done
awk -F '\t' '$1 >= 100 && $1 <= 305 && $1 ~ /^[0-9]+$/ {
  if ($1 != 100+n || $4 != "a800" || $5 != "a800" ||
      ($2 != "rfd4" && $2 != "bc2_v2")) exit 1
  n++
} END {if (n != 206) exit 1}' config/batches.tsv
groups=()
for (( start=100; start<=305; start+=GROUP_SIZE )); do groups+=("$start"); done
group_tasks=$(IFS=,; echo "${groups[*]}")
other_jobs=$(squeue -h -r -u "$USER" -o '%F' | awk -v gen="$GEN" -v col="$OLD_COL" '$1 != gen && $1 != col {n++} END {print n+0}')
(( other_jobs + 2*206 + 2*${#groups[@]} <= 500 ))
gpu_ids=$(awk -F '\t' '$1 == "generation" || $1 == "scoring" {print $2}' "$IDS" | paste -sd, -)
active_gpu=$(squeue -h -j "$gpu_ids" -o '%F' | sort -u)
if [[ -n "$active_gpu" ]]; then
  [[ "$active_gpu" == '64675209' ]]
  score_info=$(scontrol show job -o 64675209)
  grep -q 'ArrayTaskThrottle=1 ' <<< "$score_info"
fi

mkdir reports/parallel_grouped_submission.lock
scancel "$OLD_COL"
for (( check=0; check<10; check++ )); do
  [[ -z "$(squeue -h -j "$OLD_COL" -o '%i')" ]] && break
  sleep 1
done
[[ -z "$(squeue -h -j "$OLD_COL" -o '%i')" ]]
attempt=reports/parallel_grouped_attempt.tsv
printf 'generation\t%s\t3\n' "$GEN" > "$attempt"
submit_held() {
  local stage=$1 id
  shift
  id=$(sbatch --parsable --hold "$@") || return 1
  id=${id%%;*}
  [[ "$id" =~ ^[0-9]+$ ]] || return 1
  printf '%s\t%s\t3\n' "$stage" "$id" >> "$attempt"
  printf '%s\n' "$id"
}
col=$(submit_held collection --job-name=s5k_collection --array="$group_tasks%1" scripts/grouped_cpu.slurm "$ROOT" collection "$GROUP_SIZE")
score=$(submit_held scoring --array=100-305%3 scripts/scoring.slurm)
audit=$(submit_held audit --job-name=s5k_audit --array="$group_tasks%4" scripts/grouped_cpu.slurm "$ROOT" audit "$GROUP_SIZE")
for start in "${groups[@]}"; do
  gen_deps=()
  score_deps=()
  for (( task=start; task<start+GROUP_SIZE && task<=305; task++ )); do
    gen_deps+=("${GEN}_$task")
    score_deps+=("${score}_$task")
    scontrol update "JobId=${score}_$task" "Dependency=afterok:${col}_$start"
  done
  gen_dep=$(IFS=:; echo "${gen_deps[*]}")
  score_dep=$(IFS=:; echo "${score_deps[*]}")
  scontrol update "JobId=${col}_$start" "Dependency=afterany:$gen_dep"
  scontrol update "JobId=${audit}_$start" "Dependency=afterany:$score_dep,afterok:${col}_$start"
done
cat "$attempt" >> "$IDS"
printf '{"waves":[1,2,3],"batch_count":206,"cpu_group_size":%s,"cpu_groups_per_stage":%s,"generation":"%s","collection":"%s","scoring":"%s","audit":"%s","generation_concurrency":4,"scoring_concurrency":3}\n' \
  "$GROUP_SIZE" "${#groups[@]}" "$GEN" "$col" "$score" "$audit" > reports/parallel_submission.json
for job in "$audit" "$score" "$col" "$GEN"; do scontrol release "$job"; done
cat reports/parallel_submission.json
