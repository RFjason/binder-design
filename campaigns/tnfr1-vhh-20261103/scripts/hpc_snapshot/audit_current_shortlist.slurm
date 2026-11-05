#!/bin/bash
#SBATCH --job-name=tnfr1_shortlist_audit
#SBATCH --partition=64c512g
#SBATCH --exclude=node019
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:15:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/shortlist_audit_%j.out
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONNOUSERSITE=1
/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006/envs/esmfold2-esm-cu128/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/scripts/audit_current_shortlist.py
