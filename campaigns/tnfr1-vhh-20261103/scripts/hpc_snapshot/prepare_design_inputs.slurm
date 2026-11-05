#!/bin/bash
#SBATCH --job-name=tnfr1_inputs
#SBATCH --partition=64c512g
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/design_inputs_%j.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
bash /dssg/home/acct-clshzh/jasoncao1280/rfd4_deploy_20261022/scripts/run.sh python "$ROOT/scripts/prepare_design_inputs.py"
