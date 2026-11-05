"""Run an installed BC2 configuration only after an explicit CUDA backend check."""
import argparse
import glob
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sysconfig

assert os.environ.get('SLURM_JOB_ID'), 'Slurm compute nodes only'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
BC2=Path('/dssg/home/acct-clshzh/jasoncao1280/BindCraft2')
parser=argparse.ArgumentParser();parser.add_argument('config');args=parser.parse_args()
config=Path(args.config)
assert config.resolve().is_relative_to(ROOT/'inputs')
libraries=sorted(glob.glob(sysconfig.get_paths()['purelib']+'/nvidia/*/lib'))
assert libraries, 'BC2 bundled NVIDIA libraries are required'
env=os.environ.copy()
env.update(JAX_PLATFORMS='cuda',XLA_PYTHON_CLIENT_PREALLOCATE='false',PYTHONNOUSERSITE='1',
           LD_LIBRARY_PATH=':'.join(libraries+[env.get('LD_LIBRARY_PATH','')]),
           BINDCRAFT_AF2_PARAMS='/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params',
           JAX_COMPILATION_CACHE_DIR=str(ROOT/'cache/bc2_jax'),
           CUDA_CACHE_PATH=str(ROOT/'cache/bc2_gpu_required'/os.environ['SLURM_JOB_ID']))
probe="import json,jax; devices=jax.devices(); assert devices and all(d.platform=='gpu' for d in devices),devices; print(json.dumps({'backend':jax.default_backend(),'devices':[str(d) for d in devices]}))"
result=subprocess.run([str(BC2/'.venv/bin/python'),'-c',probe],env=env,text=True,capture_output=True)
report={'job_id':os.environ['SLURM_JOB_ID'],'ok':result.returncode==0,
        'config':str(config),'config_sha256':hashlib.sha256(config.read_bytes()).hexdigest(),
        'runtime_library_paths':libraries,'stdout':result.stdout,'stderr':result.stderr}
(ROOT/'reports'/('bc2_required_backend_'+os.environ['SLURM_JOB_ID']+'.json')).write_text(json.dumps(report,indent=2))
print('BC2_GPU_REQUIRED_CHECK',json.dumps(report),flush=True)
if result.returncode:raise SystemExit(result.returncode)
subprocess.run([str(BC2/'.venv/bin/bindcraft'),'design',str(config)],cwd=BC2,env=env,check=True)
print('BC2_GPU_REQUIRED_GENERATION_COMPLETE',flush=True)
