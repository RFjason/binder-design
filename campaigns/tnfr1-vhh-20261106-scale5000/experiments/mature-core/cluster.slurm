#!/bin/bash
#SBATCH --job-name=tnfr1_pose21
#SBATCH --partition=64c512g
#SBATCH --array=0
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
cd "${1:?experiment root required}"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONNOUSERSITE=1
PY=../anthropic_pdl1_pipeline_20261006/envs/esmfold2-esm-cu128/bin/python
flock -w 120 reports/pose_clusters.lock "$PY" scripts/analyze.py cluster
