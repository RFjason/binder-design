#!/bin/bash
#SBATCH --job-name=tnfr1_bc2_gpu_probe
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:05:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bc2_gpu_probe_%j.out
set -euo pipefail
export OMP_NUM_THREADS=2 PYTHONNOUSERSITE=1
/dssg/home/acct-clshzh/jasoncao1280/BindCraft2/.venv/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/scripts/bc2_gpu_backend_probe.py
