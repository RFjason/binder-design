#!/bin/bash
#SBATCH --job-name=tnfr1_bc2_v2check
#SBATCH --partition=64c512g
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bc2_v2check_%j.out
set -euo pipefail
export JAX_PLATFORMS=cpu OMP_NUM_THREADS=2 PYTHONNOUSERSITE=1
export BINDCRAFT_AF2_PARAMS=/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
cd /dssg/home/acct-clshzh/jasoncao1280/BindCraft2
./.venv/bin/python "$ROOT/scripts/prepare_bc2_v2.py"
