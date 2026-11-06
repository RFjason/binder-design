"""Run one bounded compute stage, retaining metadata for every attempt."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from s5k_common import atomic_json,read_json,new_count,progress
assert os.environ.get('SLURM_JOB_ID')
ROOT=Path(__file__).resolve().parents[1];c=read_json(ROOT/'config/campaign.json')
LEGACY=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006')
RFD=Path('/dssg/home/acct-clshzh/jasoncao1280/rfd4_deploy_20261022')
RF=Path('/dssg/home/acct-clshzh/jasoncao1280/RFantibody');BC=Path('/dssg/home/acct-clshzh/jasoncao1280/BindCraft2')
parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['generation','collection','scoring','audit']);a=parser.parse_args()
task=int(os.environ.get('S5K_TASK',os.environ.get('SLURM_ARRAY_TASK_ID','-1')));b=c['batches'][task]
assert b['task_id']==task and read_json(ROOT/'reports/prepared.json')['ok']
ident=b['batch_id'];label='scale5000_'+ident;out=Path(b['output_dir']);out.mkdir(parents=True,exist_ok=True)
route=b.get('canary_route',b['method'])
env=os.environ.copy();env.update(PYTHONNOUSERSITE='1',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='2',S5K_TASK=str(task))
env['PATH']=str(RF/'.venv/bin')+':'+env['PATH']
path=ROOT/'reports/jobs'/(ident+'_'+a.stage+'.json');old=read_json(path,{})
record=dict(batch_id=ident,method=b['method'],stage=a.stage,job_id=os.environ['SLURM_JOB_ID'],started_epoch=time.time(),hostname=os.uname().nodename,state='RUNNING',commands=[],previous_attempts=old.get('previous_attempts',[])+([old] if old else []))
def save():atomic_json(path,record)
def run(cmd,cwd=ROOT,rf_env=False):
 cmd=list(map(str,cmd))
 if rf_env:cmd=['bash','-lc','module load gcc/11.2.0 cuda/12.1.1; exec "$@"','s5k_rf',*cmd]
 record['commands'].append(cmd);save();print('S5K_COMMAND',json.dumps(cmd),flush=True)
 subprocess.run(cmd,cwd=cwd,env=env,check=True)
def budget():
 hours=c['baseline_allocated_gpu_hours']
 for p in (ROOT/'reports/jobs').glob('*.json'):
  v=read_json(p,{})
  for v in v.get('previous_attempts',[])+[v]:
   if v.get('stage') in ('generation','scoring'):
    hours+=min(2700 if v['stage']=='generation' else 1800,v.get('elapsed_seconds',max(0,time.time()-v['started_epoch'])))/3600
 ledger=read_json(ROOT/'reports/slurm_cost.json',{})
 hours=max(hours,ledger.get('cumulative_allocated_gpu_hours',0))
 record['budget_observed_gpu_hours']=hours;save()
 assert hours<c['new_gpu_stage_stop_threshold'],'S5K_BUDGET_STOP'
save()
try:
 if a.stage in ('generation','scoring'):budget()
 if a.stage=='generation':
  if task not in c['canary_tasks'].values():assert read_json(ROOT/'reports'/('route_ready_'+route+'.json'),{}).get('ok'),'Route canary not ready'
  if new_count(read_json(ROOT/'state/registry.json',{}))>=c['new_quota']:record['state']='SKIPPED_QUOTA_FULL'
  elif (out/'GENERATION_COMPLETE.json').exists():record['state']='SKIPPED_ALREADY_COMPLETE'
  else:
   if b['method']=='rfd4':
    project=RFD/'src/RFD4-Proteina';env['PROJECT_ROOT']=str(project)
    run(['bash',RFD/'scripts/run.sh','python',ROOT/'scripts/gpu_guard.py',project/'src/rfproteina/generate.py','--config-path',project/'examples/antibody_design','--config-name','vhh',f'hf_ckpt_path={RFD}/cache/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502','inference_experiment=quickstart','recycle_strategy=high/front_loaded','inputs='+b['input_cif'],'num_replicates=20','batch_size=1','num_workers=0','seed='+str(b['seed']),'compile_lengths=null','generation.dataset.sampling_transforms=null','metrics_tags=[]','out_dir='+str(out),f'hydra.run.dir={ROOT}/cache/hydra/{ident}_{os.environ["SLURM_JOB_ID"]}','prevalidate_inputs=true'])
   elif b['method']=='bc2_v2':
    if b.get('framework_contract'):run([sys.executable,ROOT/'scripts/prepare_multiframe_gpu.py',json.dumps(b)])
    run([BC/'.venv/bin/python',ROOT/'scripts/run_bc2_gpu_required.py',b['config']])
   else:
    mpnn=[RF/'.venv/bin/python',ROOT/'scripts/gpu_guard.py',RF/'scripts/proteinmpnn_interface_design.py','-checkpoint_path',RF/'weights/ProteinMPNN_v48_noise_0.2.pt','-temperature',str(b.get('temperature',.2)),'-seqs_per_struct',str(20 if b['method']=='mpnn_refine' else 4),'-loop_string','H1,H2,H3','-omit_AAs','CX','-deterministic','-debug']
    if b['method']=='rfantibody':
     hot=(ROOT/'inputs'/f'hotspots_{b["patch"]}.txt').read_text().strip()
     # CLI children inherit the explicit environment, and native models are CUDA-checked first.
     run([RF/'.venv/bin/python','-c',"import torch; assert torch.cuda.is_available(); print(torch.cuda.get_device_name())"],rf_env=True)
     run([RF/'.venv/bin/rfdiffusion','--target',ROOT/'inputs/TNFR1_mature_canonical.pdb','--framework',RF/'scripts/examples/example_inputs/h-NbBCII10.pdb','--output',out/'samples_design','--output-quiver',out/'backbones.qv','--num-designs','4','--design-loops','H1:7,H2:6,H3:10-16','--hotspots',hot,'--diffuser-t','50','--deterministic','--no-trajectory','--extra','inference.design_startnum='+str(b['seed'])],cwd=out,rf_env=True)
     run(mpnn+['-quiver',out/'backbones.qv','-outquiver',out/'sequences.qv'],cwd=out,rf_env=True)
     run([RF/'.venv/bin/qvextract',out/'sequences.qv','-o',out/'design_pdbs'],cwd=out,rf_env=True)
     run([RF/'.venv/bin/rf2','--input-quiver',out/'sequences.qv','--output-quiver',out/'rf2.qv','--hotspot-show-prop','0.0','--num-recycles','10','--seed',str(b['rf2_seed'])],cwd=out,rf_env=True)
    else:run(mpnn+['-pdbdir',b['parent_input'],'-outpdbdir',out/'design_pdbs'],cwd=out,rf_env=True)
   atomic_json(out/'GENERATION_COMPLETE.json',{'job_id':os.environ['SLURM_JOB_ID'],'batch_id':ident})
 elif a.stage=='collection':run([sys.executable,ROOT/'scripts/collect_scale5000.py'])
 elif a.stage=='scoring':
  manifest=ROOT/'inputs'/('canary_scoring_'+label+'.json')
  if read_json(manifest):
   for attempt in range(2):
    try:run(['bash',ROOT/'scripts/predict.slurm',manifest]);break
    except subprocess.CalledProcessError:
     if attempt:raise
     print('S5K_SCORE_RETRY_SAME_SEED',flush=True)
  else:record['state']='SKIPPED_NO_UNIQUE_CANDIDATES'
 else:
  run([sys.executable,ROOT/'scripts/audit_canary_scores.py',b['method'],'--batch-id',ident,'--campaign','scale5000'])
  result=read_json(ROOT/'reports'/('canary_result_'+label+'.json'))
  if task in c['canary_tasks'].values():
   assert result['pipeline_pass'] and result['scored_candidates'] and not result['failures'],'Canary requires a novel scored candidate'
   atomic_json(ROOT/'reports'/('route_ready_'+route+'.json'),{'ok':True,'canary_task':task,'job_id':os.environ['SLURM_JOB_ID'],'scored':len(result['scored_candidates'])})
 if record['state']=='RUNNING':record['state']='COMPLETED'
except Exception as exc:record.update(state='FAILED',error=repr(exc),traceback=traceback.format_exc())
finally:
 record.update(ended_epoch=time.time(),elapsed_seconds=time.time()-record['started_epoch']);save()
 if a.stage in ('collection','audit'):progress(ROOT,c)
print('S5K_STAGE_RESULT',json.dumps(record),flush=True)
if record['state']=='FAILED':raise SystemExit(1)
