"""Respect MaxSubmitJobs=500 using <=100-batch sequential waves."""
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
ROOT='/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000'
def write(name,text):(HERE/'scripts'/name).write_bytes(text.replace('REMOTE_ROOT',ROOT).encode())
write('submit.sh','''#!/bin/bash
set -euo pipefail
ROOT=REMOTE_ROOT
mkdir "$ROOT/reports/submission.lock"
IDS="$ROOT/reports/submission_ids.tsv"
: > "$IDS"
PREP=$(sbatch --parsable --hold "$ROOT/scripts/prepare.slurm")
printf 'prepare\\t%s\\t-1\\n' "$PREP" >> "$IDS"
bash "$ROOT/scripts/submit_wave.sh" 0 "$PREP"
scontrol release "$PREP"
MON=$(sbatch --parsable "$ROOT/scripts/monitor.slurm")
printf 'monitor\\t%s\\t-1\\n' "$MON" >> "$IDS"
cat "$IDS"
''')
write('submit_wave.sh','''#!/bin/bash
# Slurm management only; called on login or by the CPU controller.
set -euo pipefail
ROOT=REMOTE_ROOT
WAVE=$1
PREP=$2
case "$WAVE" in
  0) START=0; END=99 ;;
  1) START=100; END=199 ;;
  2) START=200; END=273 ;;
  3) START=274; END=305 ;;
  *) exit 2 ;;
esac
mkdir "$ROOT/reports/wave_${WAVE}.lock"
IDS="$ROOT/reports/submission_ids.tsv"
TASKS=''
while IFS=$'\\t' read -r task method ident genpart scorepart canary; do
  (( task >= START && task <= END )) || continue
  if (( WAVE > 0 )); then
    [[ -f "$ROOT/reports/route_ready_${method}.json" ]] || continue
  fi
  TASKS="${TASKS:+$TASKS,}$task"
done < "$ROOT/config/batches.tsv"
if [[ -z "$TASKS" ]]; then
  printf 'empty_wave\\t0\\t%s\\n' "$WAVE" >> "$IDS"
  exit 0
fi
GEN=$(sbatch --parsable --hold --array="$TASKS%6" "$ROOT/scripts/generation.slurm")
printf 'generation\\t%s\\t%s\\n' "$GEN" "$WAVE" >> "$IDS"
COL=$(sbatch --parsable --hold --array="$TASKS%1" "$ROOT/scripts/collection.slurm")
printf 'collection\\t%s\\t%s\\n' "$COL" "$WAVE" >> "$IDS"
SCORE=$(sbatch --parsable --hold --array="$TASKS%2" "$ROOT/scripts/scoring.slurm")
printf 'scoring\\t%s\\t%s\\n' "$SCORE" "$WAVE" >> "$IDS"
AUDIT=$(sbatch --parsable --hold --array="$TASKS%4" "$ROOT/scripts/audit.slurm")
printf 'audit\\t%s\\t%s\\n' "$AUDIT" "$WAVE" >> "$IDS"
while IFS=$'\\t' read -r task method ident genpart scorepart canary; do
  [[ ",$TASKS," == *",$task,"* ]] || continue
  DEP=''
  if (( WAVE == 0 )); then DEP="afterok:$PREP"; fi
  if (( WAVE == 0 )) && [[ "$task" != "$canary" ]]; then DEP="$DEP,afterok:${AUDIT}_$canary"; fi
  scontrol update "JobId=${GEN}_$task" "Partition=$genpart" "Dependency=$DEP"
  if (( WAVE == 0 )); then CDEP="afterok:$PREP,afterany:${GEN}_$task"; else CDEP="afterany:${GEN}_$task"; fi
  scontrol update "JobId=${COL}_$task" "Dependency=$CDEP"
  scontrol update "JobId=${SCORE}_$task" "Partition=$scorepart" "Dependency=afterok:${COL}_$task"
  scontrol update "JobId=${AUDIT}_$task" "Dependency=afterany:${SCORE}_$task,afterok:${COL}_$task"
done < "$ROOT/config/batches.tsv"
printf '{"wave":%s,"generation":"%s","collection":"%s","scoring":"%s","audit":"%s","generation_concurrency":6,"scoring_concurrency":2}\\n' "$WAVE" "$GEN" "$COL" "$SCORE" "$AUDIT" > "$ROOT/reports/wave_${WAVE}.json"
for JOB in "$AUDIT" "$SCORE" "$COL" "$GEN"; do scontrol release "$JOB"; done
cat "$ROOT/reports/wave_${WAVE}.json"
''')
print('Waves: 0..99,100..199,200..273,274..305; <=402 simultaneously submitted tasks')
