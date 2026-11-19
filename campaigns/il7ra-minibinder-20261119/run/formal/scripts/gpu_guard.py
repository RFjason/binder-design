"""Run the RFD4 entry point only after verifying CUDA inside its container."""
import os
import runpy
import sys

import torch

assert os.environ.get("SLURM_JOB_ID") and torch.cuda.is_available(), "CUDA compute node required"
print("RFD4_CUDA_REQUIRED", torch.cuda.get_device_name(), flush=True)
script = sys.argv[1]
sys.argv = sys.argv[1:]
runpy.run_path(script, run_name="__main__")
