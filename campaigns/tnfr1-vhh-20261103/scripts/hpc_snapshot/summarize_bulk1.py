"""Build an atomic cumulative snapshot from completed batch reports."""
import hashlib
import json
import os
from pathlib import Path
import time

assert os.environ.get('SLURM_JOB_ID'), 'Compute node only'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
campaign=json.loads((ROOT/'inputs/bulk1/campaign.json').read_text())
summary={'campaign':'bulk1','updated_epoch':time.time(),'batches_planned':len(campaign['batches']),
         'expected_attempts':940,'sequence_upper_bound':1600,'collection_completed':0,
         'new_unique_collected':0,'source_rejected':0,'duplicates_skipped':0,
         'scoring_audits_written':0,'scored_and_audited':0,'quality_pass':0,
         'final_submission_eligible':0,'methods':{},'stage_failures':[],'top_candidates':[],
         'completed_gpu_hours_from_stage_metadata':0.0,
         'gpu_hours_note':'Metadata sums completed generation/scoring elapsed time only; Slurm accounting is authoritative and also captures interrupted jobs.'}
all_rows=[]
for b in campaign['batches']:
    method=b['method'];label='bulk1_'+b['batch_id']
    m=summary['methods'].setdefault(method,{'unique_collected':0,'scored':0,'quality_pass':0,'max_ipsae':None})
    for stage in ('generation','collection','scoring','audit'):
        p=ROOT/'reports/bulk1_jobs'/(b['batch_id']+'_'+stage+'.json')
        if p.exists():
            try:d=json.loads(p.read_text())
            except ValueError:continue
            if d['state']=='FAILED':summary['stage_failures'].append({'batch_id':b['batch_id'],'stage':stage,'error':d.get('error')})
            for attempt in d.get('previous_attempts',[])+[d]:
                if stage in ('generation','scoring') and attempt.get('ended_epoch'):
                    summary['completed_gpu_hours_from_stage_metadata']+=attempt['elapsed_seconds']/3600
    p=ROOT/'reports'/('canary_collection_'+label+'.json')
    if p.exists():
        try:d=json.loads(p.read_text())
        except ValueError:continue
        summary['collection_completed']+=1;summary['new_unique_collected']+=len(d['candidates'])
        summary['source_rejected']+=len(d['rejected_files']);summary['duplicates_skipped']+=len(d.get('duplicates_skipped',[]))
        m['unique_collected']+=len(d['candidates'])
    p=ROOT/'reports'/('canary_result_'+label+'.json')
    if p.exists():
        try:d=json.loads(p.read_text())
        except ValueError:continue
        summary['scoring_audits_written']+=1;summary['scored_and_audited']+=len(d['scored_candidates'])
        summary['quality_pass']+=d['computational_quality_pass_count'];m['scored']+=len(d['scored_candidates'])
        m['quality_pass']+=d['computational_quality_pass_count'];all_rows.extend(d['scored_candidates'])
        if d['scored_candidates']:
            value=max(row['ipsae_min'] for row in d['scored_candidates'])
            m['max_ipsae']=value if m['max_ipsae'] is None else max(m['max_ipsae'],value)
        for f in d['failures']:summary['stage_failures'].append({'batch_id':b['batch_id'],'stage':'candidate_audit',**f})
ranked=sorted(all_rows,key=lambda r:(r['computational_quality_pass'],r['ipsae_min']),reverse=True)
for r in ranked[:50]:
    summary['top_candidates'].append({k:r.get(k) for k in ('design_id','method','stratum','parent_id','length','ipsae_min','monomer_plddt','pose_dockq','computational_quality_pass','final_rank_eligible')})
summary['distinct_scored_parent_ids']=len({r['parent_id'] for r in all_rows})
p=ROOT/'reports/bulk1_progress.json';tmp=p.with_name(p.name+'.'+os.environ['SLURM_JOB_ID']+'.tmp')
tmp.write_text(json.dumps(summary,indent=2));tmp.replace(p)
print(json.dumps({k:v for k,v in summary.items() if k not in ('top_candidates','stage_failures')}),flush=True)
