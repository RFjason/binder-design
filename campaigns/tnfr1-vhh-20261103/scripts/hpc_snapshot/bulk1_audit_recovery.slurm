#!/bin/bash
#SBATCH --job-name=tnfr1_audit_recovery
#SBATCH --partition=64c512g
#SBATCH --exclude=node019
#SBATCH --array=0-4%2
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/audit_recovery_%A_%a.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
PY=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006/envs/esmfold2-esm-cu128/bin/python
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONNOUSERSITE=1
BATCHES=(rfantibody_crd2_00 rfantibody_crd3_01 rfantibody_crd2_04 rfd4_crd2_long_04 rfd4_crd3_compact_03)
METHODS=(rfantibody rfantibody rfantibody rfd4 rfd4)
status=0
"$PY" "$ROOT/scripts/audit_canary_scores.py" "${METHODS[$SLURM_ARRAY_TASK_ID]}" --batch-id "${BATCHES[$SLURM_ARRAY_TASK_ID]}" || status=$?
"$PY" "$ROOT/scripts/summarize_bulk1.py"
exit "$status"
