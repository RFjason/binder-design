#!/bin/bash
#SBATCH --job-name=s5k_prepare
#SBATCH --partition=64c512g
#SBATCH --exclude=node019
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000/logs/prepare_%A_%a.out
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006/envs/esmfold2-esm-cu128/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000/scripts/prepare_remote.py
module load gcc/11.2.0 cuda/12.1.1
/dssg/home/acct-clshzh/jasoncao1280/RFantibody/.venv/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000/scripts/validate_parent_native.py
module unload cuda/12.1.1 gcc/11.2.0
export JAX_PLATFORMS=cpu BINDCRAFT_AF2_PARAMS=/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params
/dssg/home/acct-clshzh/jasoncao1280/BindCraft2/.venv/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000/scripts/preflight_bc2.py
