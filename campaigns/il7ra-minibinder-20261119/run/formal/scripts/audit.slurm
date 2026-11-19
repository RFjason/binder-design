#!/bin/bash
#SBATCH --job-name=il7ra5k_audit
#SBATCH --partition=64c512g
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=02:00:00
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
LEGACY=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006
cd "$SLURM_SUBMIT_DIR"
"$LEGACY/envs/esmfold2-esm-cu128/bin/python" scripts/audit_pilot.py
"$LEGACY/envs/esmfold2-esm-cu128/bin/python" scripts/cluster_scale.py
exec "$LEGACY/envs/esmfold2-esm-cu128/bin/python" scripts/select_scale.py
