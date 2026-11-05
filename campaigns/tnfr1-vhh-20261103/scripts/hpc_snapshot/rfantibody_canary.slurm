#!/bin/bash
#SBATCH --job-name=tnfr1_rfab_canary
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --output=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103/logs/rfab_canary_%j.out
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
RFAB=/dssg/home/acct-clshzh/jasoncao1280/RFantibody
OUT="$ROOT/generation/rfantibody_canary"
module load gcc/11.2.0 cuda/12.1.1
export OMP_NUM_THREADS=4 PYTHONNOUSERSITE=1
export PATH="$RFAB/.venv/bin:$PATH"
test -x "$RFAB/.venv/bin/python"
mkdir -p "$OUT"
cd "$RFAB"
HOTSPOTS=$(cat "$ROOT/inputs/hotspots_crd2.txt")
"$RFAB/.venv/bin/rfdiffusion" --target "$ROOT/inputs/TNFR1_mature_canonical.pdb" \
  --framework "$RFAB/scripts/examples/example_inputs/h-NbBCII10.pdb" \
  --output-quiver "$OUT/backbones.qv" --num-designs 1 \
  --design-loops H1:7,H2:6,H3:10-16 --hotspots "$HOTSPOTS" --diffuser-t 50 --deterministic
"$RFAB/.venv/bin/proteinmpnn" --input-quiver "$OUT/backbones.qv" \
  --output-quiver "$OUT/sequences.qv" --seqs-per-struct 2 --temperature 0.2
"$RFAB/.venv/bin/qvextract" "$OUT/sequences.qv" -o "$OUT/design_pdbs"
"$RFAB/.venv/bin/rf2" --input-quiver "$OUT/sequences.qv" \
  --output-quiver "$OUT/rf2.qv" --hotspot-show-prop 0.0 --num-recycles 10
"$RFAB/.venv/bin/qvextract" "$OUT/rf2.qv" -o "$OUT/rf2_pdbs"
echo RFANTIBODY_CANARY_OK
