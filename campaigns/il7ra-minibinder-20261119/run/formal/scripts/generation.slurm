#!/bin/bash
#SBATCH --job-name=il7ra_generate
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:25:00
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
LEGACY=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006
exec "$LEGACY/envs/esmfold2-esm-cu128/bin/python" "$SLURM_SUBMIT_DIR/scripts/run_generation.py" "$1"
