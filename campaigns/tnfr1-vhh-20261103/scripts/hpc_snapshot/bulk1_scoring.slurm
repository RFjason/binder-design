#!/bin/bash
#SBATCH --job-name=tnfr1_bulk1_score
#SBATCH --partition=a100
#SBATCH --array=0-89%2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bulk1_score_%A_%a.out
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006/envs/esmfold2-esm-cu128/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/scripts/run_bulk1.py scoring
