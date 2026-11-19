"""Run one planned pilot task using the existing installations."""
import glob
import json
import os
from pathlib import Path
import subprocess
import sys

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
HOME = Path("/dssg/home/acct-clshzh/jasoncao1280")
ROOT = Path(__file__).resolve().parents[1]
method = sys.argv[1]
batches = json.loads((ROOT / "inputs/design_batches.json").read_text())
batch = [b for b in batches if b["method"] == method][int(os.environ["SLURM_ARRAY_TASK_ID"])]
config = ROOT / "inputs" / (batch["id"] + ".json")
out = ROOT / "generation" / batch["id"]
job = os.environ["SLURM_JOB_ID"]

if method == "rfd4":
    rfd = HOME / "rfd4_deploy_20261022"
    project = rfd / "src/RFD4-Proteina"
    weights = rfd / "cache/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502"
    command = ["bash", rfd / "scripts/run.sh", "python", ROOT / "scripts/gpu_guard.py",
               project / "src/rfproteina/generate.py", f"hf_ckpt_path={weights}",
               "inference_experiment=quickstart", "recycle_strategy=high/front_loaded",
               f"inputs={config}", f"num_replicates={batch['replicates']}",
               "batch_size=1", "num_workers=0", f"seed={batch['seed']}", "compile_lengths=null",
               "metrics_tags=[]", f"out_dir={out}", f"hydra.run.dir={ROOT}/cache/hydra/{job}",
               "prevalidate_inputs=true"]
    subprocess.run(command, check=True)
    assert len(list((out / "inference_samples").rglob("*_generated.cif"))) == batch["replicates"]
elif method == "bc2":
    bc2 = HOME / "BindCraft2"
    libraries = sorted(glob.glob(str(bc2 / ".venv/lib/python*/site-packages/nvidia/*/lib")))
    assert libraries, "BindCraft2 bundled CUDA libraries required"
    env = os.environ.copy()
    env.update(JAX_PLATFORMS="cuda", XLA_PYTHON_CLIENT_PREALLOCATE="false",
               LD_LIBRARY_PATH=":".join(libraries + [env.get("LD_LIBRARY_PATH", "")]),
               BINDCRAFT_AF2_PARAMS=str(HOME / "BindCraft/params"),
               JAX_COMPILATION_CACHE_DIR=str(ROOT / "cache/bc2_jax"),
               CUDA_CACHE_PATH=str(ROOT / "cache/bc2_cuda" / job))
    probe = ("import jax; devices = jax.devices(); assert devices and all(d.platform == 'gpu' for d in devices);"
             "print('BC2_CUDA_REQUIRED', devices, flush=True)")
    subprocess.run([bc2 / ".venv/bin/python", "-c", probe], env=env, check=True)
    subprocess.run([bc2 / ".venv/bin/bindcraft", "design", config], cwd=bc2, env=env, check=True)
    assert out.is_dir()
else:
    raise ValueError(method)

(ROOT / "reports" / ("done_" + batch["id"] + ".json")).write_text(
    json.dumps({**batch, "job": job, "status": "COMPLETED"}))
