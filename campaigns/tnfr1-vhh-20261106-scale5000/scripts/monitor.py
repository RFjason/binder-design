"""Slurm CPU controller: own-job cost accounting, circuit breakers and quota stop."""
import json
import os
from pathlib import Path
import re
import subprocess
import time
from s5k_common import atomic_json,read_json,progress
assert os.environ.get('SLURM_JOB_ID')
ROOT=Path(__file__).resolve().parents[1];c=read_json(ROOT/'config/campaign.json')
ids=[];gpuids=[]
def cmd(args):return subprocess.run(args,capture_output=True,text=True,check=True).stdout
def cancel(job):
 assert any(job.startswith(i+'_') for i in ids)
 subprocess.run(['scancel',job],check=False)
events=[];start=time.monotonic()
while time.monotonic()-start<5.5*3600:
 try:
  submitted=[line.split('\t') for line in (ROOT/'reports/submission_ids.tsv').read_text().splitlines()]
  ids=[job for stage,job,wave in submitted if stage in ('generation','collection','scoring','audit')]
  gpuids=[job for stage,job,wave in submitted if stage in ('generation','scoring')]
  genids=[job for stage,job,wave in submitted if stage=='generation']
  prep=next(job for stage,job,wave in submitted if stage=='prepare')
  wave=max(int(w) for stage,job,w in submitted if int(w)>=0)
  q=cmd(['squeue','-h','-r','-u',os.environ['USER'],'-o','%i|%T|%R'])
  active=[line.split('|',2) for line in q.splitlines() if any(line.startswith(i+'_') for i in ids)]
  # JobIDRaw is a physical numeric ID; JobID preserves array parent_task IDs.
  accounting=cmd(['sacct','-X','-n','-P','-j',','.join(gpuids),'--format=JobID%100,State,ElapsedRaw,AllocTRES'])
  seconds=0;jobs=[]
  for line in accounting.splitlines():
   parts=line.split('|')
   if len(parts)<4:continue
   job,state,elapsed,tres=[p.strip() for p in parts[:4]]
   if not any(job.startswith(i+'_') for i in gpuids) or '.' in job:continue
   gpu=re.search(r'(?:^|,)gres/gpu=(\d+)(?:,|$)',tres)
   if gpu:
    sec=int(elapsed or 0)*int(gpu[1]);seconds+=sec;jobs.append({'job_id':job,'state':state,'allocated_gpu_seconds':sec})
  total=c['baseline_allocated_gpu_hours']+seconds/3600
  atomic_json(ROOT/'reports/slurm_cost.json',{'updated_epoch':time.time(),'new_allocated_gpu_hours':seconds/3600,'cumulative_allocated_gpu_hours':total,'jobs':jobs})
  p=progress(ROOT,c)
  failed_methods=set()
  for method,task in c['canary_tasks'].items():
   f=read_json(ROOT/'reports/jobs'/(c['batches'][task]['batch_id']+'_audit.json'))
   if f and f['state']=='FAILED':failed_methods.add(method)
  for method in c['canary_tasks']:
   records=[read_json(p,{}) for p in (ROOT/'reports/jobs').glob('*_generation.json')]
   records=sorted([v for v in records if v.get('method')==method and v.get('state') in ('FAILED','COMPLETED')],key=lambda v:v.get('ended_epoch',0))
   if len(records)>=3 and all(v['state']=='FAILED' for v in records[-3:]):failed_methods.add(method)
  for job,state,reason in active:
   if state!='PENDING':continue
   arr,_,task=job.partition('_')
   if not task.isdigit():continue
   task=int(task)
   isgen=arr in genids
   stop=(isgen and (p['new_reserved']>=c['new_quota'] or c['batches'][task].get('canary_route',c['batches'][task]['method']) in failed_methods))
   stop=stop or (total>=c['new_gpu_stage_stop_threshold'] and arr in gpuids)
   stop=stop or reason.startswith('DependencyNeverSatisfied')
   if stop:
    cancel(job);events.append({'epoch':time.time(),'cancelled_pending':job,'reason':reason,'quota':p['new_reserved'],'gpu_hours':total})
  atomic_json(ROOT/'reports/monitor.json',{'job_id':os.environ['SLURM_JOB_ID'],'updated_epoch':time.time(),'active_tasks':len(active),'failed_methods':sorted(failed_methods),'events':events})
  if not active:
   if wave<3 and p['new_reserved']<c['new_quota'] and total<c['new_gpu_stage_stop_threshold']:
    # Completed Slurm jobs can age out. Later waves use the recorded native
    # preflight and route-ready files instead of old scheduler dependencies.
    result=cmd(['bash',str(ROOT/'scripts/submit_wave.sh'),str(wave+1),prep])
    print('S5K_NEXT_WAVE',result,flush=True)
    continue
   atomic_json(ROOT/'reports/campaign_finished.json',{'updated_epoch':time.time(),'progress':p,'target_reached':p['cumulative_scored']>=c['target'],'gpu_hours':total})
   print('S5K_FINISHED',json.dumps(p),flush=True);break
 except Exception as exc:
  events.append({'epoch':time.time(),'monitor_error':repr(exc)})
 time.sleep(60)
else:
 job=cmd(['sbatch','--parsable',str(ROOT/'scripts/monitor.slurm')]).strip()
 atomic_json(ROOT/'reports/monitor_continuation.json',{'from_job':os.environ['SLURM_JOB_ID'],'next_job':job,'epoch':time.time()})
 print('S5K_MONITOR_CONTINUED',job,flush=True)
