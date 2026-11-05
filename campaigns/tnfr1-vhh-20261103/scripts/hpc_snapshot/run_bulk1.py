"""Execute one finite batch stage on its assigned Slurm compute node."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

assert os.environ.get('SLURM_JOB_ID'), 'All stages require a Slurm compute node'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
LEGACY=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006')
RFD4=Path('/dssg/home/acct-clshzh/jasoncao1280/rfd4_deploy_20261022')
RFAB=Path('/dssg/home/acct-clshzh/jasoncao1280/RFantibody')
BC2=Path('/dssg/home/acct-clshzh/jasoncao1280/BindCraft2')
campaign=json.loads((ROOT/'inputs/bulk1/campaign.json').read_text())
preflight=json.loads((ROOT/'reports/bulk1_preflight.json').read_text())
assert preflight['ok']
assert hashlib.sha256((ROOT/'inputs/bulk1/campaign.json').read_bytes()).hexdigest()==preflight['hashes'][str(ROOT/'inputs/bulk1/campaign.json')]
parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['generation','collection','scoring','audit'])
args=parser.parse_args()
task=int(os.environ['SLURM_ARRAY_TASK_ID'])
b=campaign['batches'][task]
assert b['task_id']==task
ident=b['batch_id'];label='bulk1_'+ident;out=Path(b['output_dir'])
env=os.environ.copy();env.update(PYTHONNOUSERSITE='1',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='4')
env['PATH']=str(RFAB/'.venv/bin')+':'+env['PATH']
meta_dir=ROOT/'reports/bulk1_jobs';meta_dir.mkdir(parents=True,exist_ok=True)
meta_file=meta_dir/(ident+'_'+args.stage+'.json')
previous=[]
if meta_file.exists():
    old=json.loads(meta_file.read_text())
    previous=old.pop('previous_attempts',[])+[old]
meta={'batch':b,'stage':args.stage,'job_id':os.environ['SLURM_JOB_ID'],
      'array_job_id':os.environ.get('SLURM_ARRAY_JOB_ID'),'task_id':task,
      'hostname':os.uname().nodename,'started_epoch':time.time(),'state':'RUNNING','commands':[],
      'previous_attempts':previous}
def save():
    tmp=meta_file.with_suffix('.tmp');tmp.write_text(json.dumps(meta,indent=2));tmp.replace(meta_file)
def run(cmd,cwd=ROOT):
    cmd=list(map(str,cmd));meta['commands'].append(cmd);save()
    print('BULK_COMMAND',json.dumps(cmd),flush=True)
    subprocess.run(cmd,check=True,cwd=cwd,env=env)
save()
try:
    if args.stage in ('generation','scoring'):
        consumed=0.0
        for p in meta_dir.glob('*.json'):
            record=json.loads(p.read_text())
            for attempt in record.get('previous_attempts',[])+[record]:
                if attempt['stage'] not in ('generation','scoring'):continue
                limit=7200 if attempt['stage']=='generation' else 3600
                consumed+=min(limit,attempt.get('elapsed_seconds',max(0,time.time()-attempt['started_epoch'])))/3600
        # Reserve all eight active workers' maximum remaining walltime plus startup overhead.
        threshold=campaign['gpu_hours_budget']-(6*2+2*1+2)
        meta['budget_check']={'observed_gpu_hours':consumed,'new_batch_stop_threshold':threshold,
                              'campaign_gpu_hours_budget':campaign['gpu_hours_budget']}
        assert consumed < threshold, 'BULK1_RESOURCE_BUDGET_STOP'
    if args.stage=='generation':
        out.mkdir(parents=True,exist_ok=True)
        if (out/'GENERATION_COMPLETE.json').exists():
            meta['state']='SKIPPED_ALREADY_COMPLETE'
        elif b['method']=='rfd4':
            project=RFD4/'src/RFD4-Proteina';env['PROJECT_ROOT']=str(project)
            run(['bash',RFD4/'scripts/run.sh','python',project/'src/rfproteina/generate.py',
                 '--config-path',project/'examples/antibody_design','--config-name','vhh',
                 f'hf_ckpt_path={RFD4}/cache/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502',
                 'inference_experiment=quickstart','recycle_strategy=high/front_loaded',
                 'inputs='+b['input_cif'],'num_replicates=20','batch_size=1','num_workers=0',
                 'seed='+str(b['seed']),'compile_lengths=null','generation.dataset.sampling_transforms=null',
                 'metrics_tags=[]','out_dir='+str(out),
                 f'hydra.run.dir={ROOT}/cache/hydra/bulk1_{ident}_{os.environ["SLURM_JOB_ID"]}',
                 'prevalidate_inputs=true'])
        elif b['method']=='rfantibody':
            hot=(ROOT/'inputs'/f'hotspots_{b["patch"]}.txt').read_text().strip()
            run([RFAB/'.venv/bin/rfdiffusion','--target',ROOT/'inputs/TNFR1_mature_canonical.pdb',
                 '--framework',RFAB/'scripts/examples/example_inputs/h-NbBCII10.pdb',
                 '--output',out/'samples_design','--output-quiver',out/'backbones.qv',
                 '--num-designs','10','--design-loops','H1:7,H2:6,H3:10-16',
                 '--hotspots',hot,'--diffuser-t','50','--deterministic','--no-trajectory',
                 '--extra','inference.design_startnum='+str(b['seed'])],cwd=RFAB)
            run([RFAB/'.venv/bin/proteinmpnn','--input-quiver',out/'backbones.qv',
                 '--output-quiver',out/'sequences.qv','--seqs-per-struct','4',
                 '--temperature','0.2','--deterministic'],cwd=RFAB)
            run([RFAB/'.venv/bin/qvextract',out/'sequences.qv','-o',out/'design_pdbs'],cwd=RFAB)
            run([RFAB/'.venv/bin/rf2','--input-quiver',out/'sequences.qv',
                 '--output-quiver',out/'rf2.qv','--hotspot-show-prop','0.0',
                 '--num-recycles','10','--seed',str(b['rf2_seed'])],cwd=RFAB)
            run([RFAB/'.venv/bin/qvextract',out/'rf2.qv','-o',out/'rf2_pdbs'],cwd=RFAB)
        else:
            env.update(XLA_PYTHON_CLIENT_PREALLOCATE='false',
                       JAX_COMPILATION_CACHE_DIR=str(ROOT/'cache/bc2_jax'),
                       CUDA_CACHE_PATH=str(ROOT/'cache/bulk1_cuda'/ident),
                       BINDCRAFT_AF2_PARAMS='/dssg/home/acct-clshzh/jasoncao1280/BindCraft/params')
            run([BC2/'.venv/bin/bindcraft','design',b['config']],cwd=BC2)
        (out/'GENERATION_COMPLETE.json').write_text(json.dumps({'batch_id':ident,'job_id':os.environ['SLURM_JOB_ID']}))
    elif args.stage=='collection':
        run([sys.executable,ROOT/'scripts/collect_canary.py',b['method'],'--batch-id',ident])
    elif args.stage=='scoring':
        manifest=ROOT/'inputs'/('canary_scoring_'+label+'.json')
        if json.loads(manifest.read_text()):
            run(['bash',ROOT/'scripts/predict.slurm',manifest])
        else:
            meta['state']='SKIPPED_NO_UNIQUE_VALID_CANDIDATES'
    else:
        run([sys.executable,ROOT/'scripts/audit_canary_scores.py',b['method'],'--batch-id',ident])
    if meta['state']=='RUNNING':meta['state']='COMPLETED'
except Exception as exc:
    meta.update(state='FAILED',error=repr(exc),traceback=traceback.format_exc())
finally:
    meta['ended_epoch']=time.time();meta['elapsed_seconds']=meta['ended_epoch']-meta['started_epoch'];save()
    if args.stage=='audit':
        subprocess.run([sys.executable,str(ROOT/'scripts/summarize_bulk1.py')],cwd=ROOT,env=env)
print('BULK_STAGE_RESULT',json.dumps(meta),flush=True)
if meta['state']=='FAILED':raise SystemExit(1)
