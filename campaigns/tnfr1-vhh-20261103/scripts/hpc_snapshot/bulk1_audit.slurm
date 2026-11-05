#!/bin/bash
#SBATCH --job-name=tnfr1_bulk1_audit
#SBATCH --partition=64c512g
#SBATCH --exclude=node019
#SBATCH --array=0-89%4
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bulk1_audit_%A_%a.out
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006/envs/esmfold2-esm-cu128/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/scripts/run_bulk1.py audit
