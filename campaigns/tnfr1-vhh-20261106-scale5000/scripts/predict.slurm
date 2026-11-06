#!/bin/bash
#SBATCH --job-name=tnfr1_esm2
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000/logs/predict_%j.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000
LEGACY=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 OMP_NUM_THREADS=4
export HF_HOME="$LEGACY/cache/huggingface"
export ESMCFOLD_CCD_PATH="$LEGACY/weights/ESMFold2/ccd.pkl"
export CUDA_CACHE_PATH="$ROOT/cache/cuda" TRITON_CACHE_DIR="$ROOT/cache/triton" TORCHINDUCTOR_CACHE_DIR="$ROOT/cache/inductor"
mkdir -p "$CUDA_CACHE_PATH" "$TRITON_CACHE_DIR" "$TORCHINDUCTOR_CACHE_DIR"
cd "$ROOT"
"$LEGACY/envs/esmfold2-esm-cu128/bin/python" "$ROOT/scripts/predict_esmfold2.py" "$1"
