"""Finite GPU-verified BC2 recovery with local dependencies and cumulative budget control."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

assert os.environ.get('SLURM_JOB_ID'), 'Compute nodes only'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
BC2=Path('/dssg/home/acct-clshzh/jasoncao1280/BindCraft2')
campaign=json.loads((ROOT/'inputs/bc2_recovery/campaign.json').read_text())
parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['generation','collection','scoring','audit']);args=parser.parse_args()
b=campaign['batches'][int(os.environ['SLURM_ARRAY_TASK_ID'])];ident=b['batch_id'];label='bc2_recovery_'+ident
folder=ROOT/'reports/bc2_recovery_jobs';folder.mkdir(exist_ok=True)
path=folder/(ident+'_'+args.stage+'.json')
record=dict(batch_id=ident,stage=args.stage,job_id=os.environ['SLURM_JOB_ID'],started_epoch=time.time(),state='RUNNING')
def save():
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2));tmp.replace(path)
def run(command):
    print('BC2_RECOVERY_COMMAND',json.dumps(list(map(str,command))),flush=True)
    subprocess.run(list(map(str,command)),check=True,cwd=ROOT)
def summarize():
    result=dict(campaign='bc2_recovery',updated_epoch=time.time(),planned_trajectories=40,
                collected=0,scored=0,quality_pass=0,generated_batches=0,stage_failures=[],
                baseline_gpu_hours_ceiling=68.0,completed_gpu_hours=0.0)
    for batch in campaign['batches']:
        prefix='bc2_recovery_'+batch['batch_id']
        for stage in ('generation','collection','scoring','audit'):
            f=folder/(batch['batch_id']+'_'+stage+'.json')
            if f.exists():
                v=json.loads(f.read_text())
                if stage in ('generation','scoring') and 'ended_epoch' in v:result['completed_gpu_hours']+=v['elapsed_seconds']/3600
                if stage=='generation' and v['state']=='COMPLETED':result['generated_batches']+=1
                if v['state']=='FAILED':result['stage_failures'].append(v)
        f=ROOT/'reports'/('canary_collection_'+prefix+'.json')
        if f.exists():
            try:result['collected']+=len(json.loads(f.read_text())['candidates'])
            except ValueError:pass
        f=ROOT/'reports'/('canary_result_'+prefix+'.json')
        if f.exists():
            try:
                v=json.loads(f.read_text());result['scored']+=len(v['scored_candidates']);result['quality_pass']+=v['computational_quality_pass_count']
            except ValueError:pass
    target=ROOT/'reports/bc2_recovery_progress.json';temp=target.with_name(target.name+'.'+os.environ['SLURM_JOB_ID']+'.tmp')
    temp.write_text(json.dumps(result,indent=2));temp.replace(target)
save()
try:
    if args.stage in ('generation','scoring'):
        hours=campaign['baseline_allocated_gpu_hours_ceiling']
        for f in folder.glob('*.json'):
            v=json.loads(f.read_text())
            if v['stage'] not in ('generation','scoring'):continue
            limit=1800 if v['stage']=='generation' else 1200
            hours+=min(limit,v.get('elapsed_seconds',time.time()-v['started_epoch']))/3600
        record['combined_gpu_hours_observed_ceiling']=hours
        assert hours < campaign['new_gpu_stage_stop_threshold'],'BC2_RECOVERY_BUDGET_STOP'
    if args.stage=='generation':
        run([BC2/'.venv/bin/python',ROOT/'scripts/run_bc2_gpu_required.py',b['config']])
    elif args.stage=='collection':
        run([sys.executable,ROOT/'scripts/collect_canary.py','bc2_v2','--batch-id',ident,'--campaign','bc2_recovery'])
    elif args.stage=='scoring':
        f=ROOT/'inputs'/('canary_scoring_'+label+'.json')
        if json.loads(f.read_text()):run(['bash',ROOT/'scripts/predict.slurm',f])
        else:record['state']='SKIPPED_NO_VALID_UNIQUE_CANDIDATE'
    else:
        run([sys.executable,ROOT/'scripts/audit_canary_scores.py','bc2_v2','--batch-id',ident,'--campaign','bc2_recovery'])
    if record['state']=='RUNNING':record['state']='COMPLETED'
except Exception as exc:
    record.update(state='FAILED',error=repr(exc),traceback=traceback.format_exc())
finally:
    record.update(ended_epoch=time.time(),elapsed_seconds=time.time()-record['started_epoch']);save()
    if args.stage in ('collection','audit'):summarize()
print('BC2_RECOVERY_STAGE',json.dumps(record),flush=True)
if record['state']=='FAILED':raise SystemExit(1)
