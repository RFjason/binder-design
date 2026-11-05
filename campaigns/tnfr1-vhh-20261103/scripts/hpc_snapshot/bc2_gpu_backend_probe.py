"""Require BC2's own CUDA libraries and fail before design if CUDA is unavailable."""
import glob
import json
import os
from pathlib import Path
import subprocess
import sysconfig

assert os.environ.get('SLURM_JOB_ID'), 'Slurm compute nodes only'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
BC2=Path('/dssg/home/acct-clshzh/jasoncao1280/BindCraft2')
env=os.environ.copy()
libraries=sorted(glob.glob(sysconfig.get_paths()['purelib']+'/nvidia/*/lib'))
assert libraries, 'BC2 bundled NVIDIA runtime libraries not found'
env['LD_LIBRARY_PATH']=':'.join(libraries+[env.get('LD_LIBRARY_PATH','')])
env.update(JAX_PLATFORMS='cuda',XLA_PYTHON_CLIENT_PREALLOCATE='false',PYTHONNOUSERSITE='1')
probe="import json,jax; devices=jax.devices(); assert devices and all(d.platform=='gpu' for d in devices),devices; print(json.dumps({'backend':jax.default_backend(),'devices':[str(d) for d in devices]}))"
result=subprocess.run([str(BC2/'.venv/bin/python'),'-c',probe],env=env,text=True,capture_output=True)
report={'job_id':os.environ['SLURM_JOB_ID'],'ok':result.returncode==0,
        'runtime_library_paths':libraries,'jax_platforms':'cuda','stdout':result.stdout,'stderr':result.stderr}
(ROOT/'reports/bc2_gpu_backend_probe.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report),flush=True)
raise SystemExit(result.returncode)
