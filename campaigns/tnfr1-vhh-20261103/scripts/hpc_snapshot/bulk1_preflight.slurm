#!/bin/bash
#SBATCH --job-name=tnfr1_bulk1_preflight
#SBATCH --partition=64c512g
#SBATCH --exclude=node019
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bulk1_preflight_%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export JAX_PLATFORMS=cpu
export BINDCRAFT_AF2_PARAMS=/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params
/dssg/home/acct-clshzh/jasoncao1280/BindCraft2/.venv/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/scripts/preflight_bulk1.py
