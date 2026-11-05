#!/bin/bash
# Login-node Slurm management only.
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
mkdir "$ROOT/reports/bc2_recovery_submission.lock"
G=$(sbatch --parsable --hold --array=0-39%4 --job-name=tnfr1_bc2_gpu_gen "$ROOT/scripts/bc2_recovery_gpu.slurm" generation)
C=$(sbatch --parsable --hold --array=0-39%4 "$ROOT/scripts/bc2_recovery_cpu.slurm" collection)
S=$(sbatch --parsable --hold --array=0-39%2 --time=00:20:00 --job-name=tnfr1_bc2_score --kill-on-invalid-dep=yes "$ROOT/scripts/bc2_recovery_gpu.slurm" scoring)
A=$(sbatch --parsable --hold --array=0-39%4 "$ROOT/scripts/bc2_recovery_cpu.slurm" audit)
printf 'generation\t%s\ncollection\t%s\nscoring\t%s\naudit\t%s\n' "$G" "$C" "$S" "$A" > "$ROOT/reports/bc2_recovery_submission_ids.tsv"
for ((i=0;i<40;i++)); do
    genpart=a100;scorepart=a800
    if ((i%2)); then genpart=a800;scorepart=a100; fi
    scontrol update "JobId=${G}_$i" "Partition=$genpart"
    scontrol update "JobId=${C}_$i" "Dependency=afterany:${G}_$i"
    scontrol update "JobId=${S}_$i" "Partition=$scorepart" "Dependency=afterok:${C}_$i"
    scontrol update "JobId=${A}_$i" "Dependency=afterany:${S}_$i"
done
printf '{"generation":"%s","collection":"%s","scoring":"%s","audit":"%s","generation_concurrency":4,"scoring_concurrency":2}\n' "$G" "$C" "$S" "$A" > "$ROOT/reports/bc2_recovery_submission.json"
scontrol release "$A"
scontrol release "$S"
scontrol release "$C"
scontrol release "$G"
cat "$ROOT/reports/bc2_recovery_submission.json"
