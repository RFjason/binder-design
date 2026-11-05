#!/bin/bash
#SBATCH --job-name=tnfr1_bc2_canary
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bc2_canary_%j.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$ROOT/cache/bc2_jax"
export CUDA_CACHE_PATH="$ROOT/cache/bc2_cuda"
export BINDCRAFT_AF2_PARAMS=/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params
mkdir -p "$JAX_COMPILATION_CACHE_DIR" "$CUDA_CACHE_PATH"
cd /dssg/home/acct-clshzh/jasoncao1280/BindCraft2
./.venv/bin/bindcraft design "$ROOT/inputs/bc2_canary_crd2.json"
