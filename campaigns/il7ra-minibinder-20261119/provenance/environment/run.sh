#!/bin/bash
set -euo pipefail
if [[ -z ${SLURM_JOB_ID:-} ]]; then
  echo 'Run only inside a Slurm job. Submit scripts/generate.slurm with sbatch.' >&2
  exit 2
fi
ROOT=/dssg/home/acct-clshzh/jasoncao1280/rfd4_deploy_20261022
IMAGE=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006/containers/ubuntu22.04.sif
export ROOT
export HF_HOME="$ROOT/cache/huggingface"
export HF_HUB_DISABLE_XET=1
export HF_HUB_OFFLINE=1
export DATA_PATH="$ROOT/data_path"
export LD_LIBRARY_PATH="$ROOT/runtime-libs:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MPLCONFIGDIR="$ROOT/cache/matplotlib"
export TORCHINDUCTOR_CACHE_DIR="$ROOT/cache/inductor"
export TRITON_CACHE_DIR="$ROOT/cache/triton"
export TORCHDYNAMO_DISABLE=1
export HYDRA_FULL_ERROR=1
export CLUSTER=local
cd "$ROOT/src/RFD4-Proteina"
entry="$1"
shift
if [[ "$entry" == rfd4-generate ]]; then
  set -- "hf_ckpt_path=$ROOT/cache/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502" "$@"
fi
exec apptainer exec --nv --env "LD_LIBRARY_PATH=$ROOT/runtime-libs:/.singularity.d/libs" --bind /dssg:/dssg "$IMAGE" "$ROOT/venv/bin/$entry" "$@"