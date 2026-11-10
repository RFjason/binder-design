#!/bin/bash
#SBATCH --job-name=tnfr1_mature21
#SBATCH --partition=a800
#SBATCH --array=0
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
set -euo pipefail
ROOT=${1:?experiment root required}
LEGACY=$(dirname "$ROOT")/anthropic_pdl1_pipeline_20261006
PY="$LEGACY/envs/esmfold2-esm-cu128/bin/python"
cd "$ROOT"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 HF_HOME="$LEGACY/cache/huggingface"
export ESMCFOLD_CCD_PATH="$LEGACY/weights/ESMFold2/ccd.pkl"
export CUDA_CACHE_PATH="$ROOT/cache/cuda" TRITON_CACHE_DIR="$ROOT/cache/triton" TORCHINDUCTOR_CACHE_DIR="$ROOT/cache/inductor"
mkdir -p "$CUDA_CACHE_PATH" "$TRITON_CACHE_DIR" "$TORCHINDUCTOR_CACHE_DIR"
"$PY" scripts/gpu_budget.py
# If the CPU partition is still down, produce the historical clustering here too.
flock -w 120 reports/pose_clusters.lock "$PY" scripts/analyze.py cluster || echo 'Historical pose cluster incomplete'
status=0
timeout --kill-after=30s 98m "$PY" scripts/predict_esmfold2.py inputs/manifest.json || status=$?
audit_status=0
"$PY" scripts/analyze.py audit || audit_status=$?
if (( status != 0 || audit_status != 0 )); then
    echo "Partial batch retained: prediction=$status audit=$audit_status"
    exit 1
fi
