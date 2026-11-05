#!/bin/bash
# Slurm management only; no computation on the login node.
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
while IFS=$'\t' read -r task method ident genpart scorepart; do
    case "$method" in
        rfd4) limit=00:30:00 ;;
        rfantibody) limit=00:45:00 ;;
        bc2_v2) limit=00:20:00 ;;
    esac
    scontrol update "JobId=64567743_$task" "TimeLimit=$limit"
done < "$ROOT/inputs/bulk1/batches.tsv"
