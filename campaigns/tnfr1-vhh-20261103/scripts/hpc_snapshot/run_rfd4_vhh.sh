#!/bin/bash
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
RFD4=/dssg/home/acct-clshzh/jasoncao1280/rfd4_deploy_20261022
export PROJECT_ROOT="$RFD4/src/RFD4-Proteina"
# The shared rfd4-generate launcher prepends a positional override, which breaks
# Hydra when --config-path follows it. Use its unchanged Python runtime instead.
bash "$RFD4/scripts/run.sh" python "$PROJECT_ROOT/src/rfproteina/generate.py" \
  --config-path "$PROJECT_ROOT/examples/antibody_design" --config-name vhh "$@" \
  "hf_ckpt_path=$RFD4/cache/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502" \
  inference_experiment=quickstart recycle_strategy=high/front_loaded \
  "inputs=$ROOT/inputs/rfd4/tnfr1_crd2_compact.cif" \
  num_replicates=1 batch_size=1 num_workers=0 seed=101 compile_lengths=null \
  generation.dataset.sampling_transforms=null metrics_tags=[] \
  "out_dir=$ROOT/generation/rfd4_canary" \
  "hydra.run.dir=$ROOT/cache/hydra/rfd4_$SLURM_JOB_ID" prevalidate_inputs=true
