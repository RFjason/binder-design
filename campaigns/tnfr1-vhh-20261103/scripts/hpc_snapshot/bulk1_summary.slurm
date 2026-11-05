#!/bin/bash
#SBATCH --job-name=tnfr1_bulk1_summary
#SBATCH --partition=64c512g
#SBATCH --exclude=node019
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/bulk1_summary_%j.out
set -euo pipefail
/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006/envs/esmfold2-esm-cu128/bin/python /dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/scripts/summarize_bulk1.py
