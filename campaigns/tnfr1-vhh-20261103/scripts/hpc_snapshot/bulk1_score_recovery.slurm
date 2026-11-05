#!/bin/bash
#SBATCH --job-name=tnfr1_score_recovery
#SBATCH --partition=a100
#SBATCH --array=0-2%2
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:45:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/score_recovery_%A_%a.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
BATCHES=(rfantibody_crd2_04 rfd4_crd2_long_04 rfd4_crd3_compact_03)
batch=${BATCHES[$SLURM_ARRAY_TASK_ID]}
bash "$ROOT/scripts/predict.slurm" "$ROOT/inputs/canary_scoring_bulk1_$batch.json"
