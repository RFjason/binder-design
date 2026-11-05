#!/bin/bash
#SBATCH --job-name=tnfr1_score_pool
#SBATCH --partition=a100
#SBATCH --array=0-3%2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/score_pool_%A_%a.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
METHODS=(rfd4 rfantibody bc2_v2 rfantibody_crd3)
METHOD=${METHODS[$SLURM_ARRAY_TASK_ID]}
# One array enforces the two-GPU scoring cap; task dependencies are independent.
bash "$ROOT/scripts/predict.slurm" "$ROOT/inputs/canary_scoring_$METHOD.json"
