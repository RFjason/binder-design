"""Write Linux submission scripts and explicit dependencies without model execution."""
from pathlib import Path
ROOT='/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000'
LEG='/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006'
RF='/dssg/home/acct-clshzh/jasoncao1280/RFantibody'
BC='/dssg/home/acct-clshzh/jasoncao1280/BindCraft2'
HERE=Path(__file__).resolve().parents[1]
def write(name,txt):(HERE/'scripts'/name).write_bytes(txt.encode())
def slurm(name,body,gpu=False,array=None,limit='00:20:00',concurrency=None):
 s=f'#!/bin/bash\n#SBATCH --job-name=s5k_{name}\n#SBATCH --partition={"a800" if gpu else "64c512g"}\n'
 s+=('#SBATCH --gres=gpu:1\n' if gpu else '#SBATCH --exclude=node019\n')
 s+=f'#SBATCH --cpus-per-task={4 if gpu else 2}\n#SBATCH --mem={64 if gpu else 16}G\n#SBATCH --time={limit}\n#SBATCH --output={ROOT}/logs/{name}_%A_%a.out\n'
 if array:s+=f'#SBATCH --array={array}%{concurrency}\n'
 s+='set -euo pipefail\nexport PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2\n'
 write(name+'.slurm',s+body+'\n')
py=LEG+'/envs/esmfold2-esm-cu128/bin/python'
slurm('prepare',f'{py} {ROOT}/scripts/prepare_remote.py\nmodule load gcc/11.2.0 cuda/12.1.1\n{RF}/.venv/bin/python {ROOT}/scripts/validate_parent_native.py\nmodule unload cuda/12.1.1 gcc/11.2.0\nexport JAX_PLATFORMS=cpu BINDCRAFT_AF2_PARAMS=/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params\n{BC}/.venv/bin/python {ROOT}/scripts/preflight_bc2.py')
for stage,arr,conc,gpu,limit in [('generation','0-273',6,True,'00:45:00'),('reserve','274-305',6,True,'00:45:00'),('collection','0-305',1,False,'00:20:00'),('scoring','0-305',2,True,'00:30:00'),('audit','0-305',4,False,'00:20:00')]:
 slurm(stage,f'{py} {ROOT}/scripts/run_scale5000.py '+('generation' if stage=='reserve' else stage),gpu,arr,limit,conc)
slurm('monitor',f'{py} {ROOT}/scripts/monitor.py',limit='06:00:00')
import runpy
runpy.run_path(str(HERE/'scripts/write_wave_submission_local.py'),run_name='__main__')
