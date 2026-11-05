#!/bin/bash
#SBATCH --job-name=tnfr1_launchcheck
#SBATCH --partition=64c512g
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/launchcheck_%j.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
bash "$ROOT/scripts/run_rfd4_vhh.sh" --cfg job --resolve > "$ROOT/reports/rfd4_launch_config_v2.yaml"
echo RFD4_LAUNCH_CONFIG_OK
