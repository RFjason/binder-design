#!/bin/bash
#SBATCH --job-name=il7ra_score
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:40:00
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
ROOT="$SLURM_SUBMIT_DIR"
LEGACY=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 HF_HOME="$LEGACY/cache/huggingface"
export ESMCFOLD_CCD_PATH="$LEGACY/weights/ESMFold2/ccd.pkl"
export CUDA_CACHE_PATH="$ROOT/cache/cuda" TRITON_CACHE_DIR="$ROOT/cache/triton"
export TORCHINDUCTOR_CACHE_DIR="$ROOT/cache/inductor"
manifest="$1"
if [[ "$manifest" == *%a* ]]; then
    manifest="${manifest//%a/$SLURM_ARRAY_TASK_ID}"
fi
exec "$LEGACY/envs/esmfold2-esm-cu128/bin/python" "$ROOT/scripts/predict_esmfold2.py" "$ROOT/$manifest"
