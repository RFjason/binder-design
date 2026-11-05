#!/bin/bash
#SBATCH --job-name=tnfr1_bc2_inputcheck
#SBATCH --partition=64c512g
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bc2_inputcheck_%j.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
export JAX_PLATFORMS=cpu OMP_NUM_THREADS=2 PYTHONNOUSERSITE=1
export BINDCRAFT_AF2_PARAMS=/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params
cd /dssg/home/acct-clshzh/jasoncao1280/BindCraft2
./.venv/bin/python "$ROOT/scripts/preflight_bc2.py"
