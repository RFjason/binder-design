#!/bin/bash
#SBATCH --job-name=il7ra_epitopes
#SBATCH --partition=64c512g
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
LEGACY=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006
exec "$LEGACY/envs/esmfold2-esm-cu128/bin/python" "$SLURM_SUBMIT_DIR/scripts/prepare_design.py" "$@"
