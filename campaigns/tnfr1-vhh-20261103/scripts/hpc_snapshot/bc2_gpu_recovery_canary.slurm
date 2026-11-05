#!/bin/bash
#SBATCH --job-name=tnfr1_bc2_gpu_canary
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bc2_gpu_recovery_%j.out
set -euo pipefail
export OMP_NUM_THREADS=4 PYTHONNOUSERSITE=1
# BC2 supplies its own NVIDIA runtime libraries; do not load another CUDA module.
/dssg/home/acct-clshzh/jasoncao1280/BindCraft2/.venv/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/scripts/run_bc2_gpu_required.py /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/inputs/bc2_gpu_recovery_canary.json
