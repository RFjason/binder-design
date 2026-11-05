#!/bin/bash
#SBATCH --job-name=tnfr1_bc2_v2
#SBATCH --partition=a100
#SBATCH --array=0-1%2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bc2_v2_%A_%a.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
PATCHES=(crd2 crd3)
PATCH=${PATCHES[$SLURM_ARRAY_TASK_ID]}
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$ROOT/cache/bc2_jax"
export CUDA_CACHE_PATH="$ROOT/cache/bc2_cuda"
export BINDCRAFT_AF2_PARAMS=/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params
cd /dssg/home/acct-clshzh/jasoncao1280/BindCraft2
./.venv/bin/bindcraft design "$ROOT/inputs/bc2_canary_v2_$PATCH.json"
