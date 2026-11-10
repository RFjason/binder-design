"""Export complete audit rows from downloaded reports, without rerunning inference."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

parser=argparse.ArgumentParser();parser.add_argument('reports',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args()
HERE=Path(__file__).resolve().parents[1];campaign=json.loads((HERE/'config/campaign.json').read_text());baseline=json.loads((HERE/'inputs/baseline_registry.json').read_text())
args.output.mkdir(parents=True,exist_ok=True)
def write(name,data): (args.output/name).write_bytes((json.dumps(data,indent=2,ensure_ascii=False)+'\n').encode())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
rows={};seqs=set();sources={};errors=[]
for p in sorted(args.reports.glob('canary_result_scale5000_*.json')):
 data=json.loads(p.read_text());sources[p.name]=digest(p);errors.extend(data.get('failures',[]))
 for row in data['scored_candidates']:
  seq=row['sequence'];key=hashlib.sha256(seq.encode()).hexdigest();cid=row['design_id']
  assert cid not in rows and key not in seqs and key not in baseline
  assert cid.endswith(key[:12]) and row['length']==len(seq) and 100<=len(seq)<=150 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY')
  assert row['scoring_protocol_verified'] and row['final_rank_eligible'] is False
  assert row['computational_quality_pass']==(row['ipsae_min']>=.6 and row['monomer_plddt']>=70 and row['pose_dockq'] is not None and row['pose_dockq']>=.23 and row['monomer_disulfide_geometry_ok'])
  row['source_audit_report']=p.name;row['framework_id']=row.get('framework_id','h-NbBCII10')
  for f in ('source_file','reference_pdb'):
   if row.get(f,'').startswith(campaign['remote_root']+'/'):row[f]=row[f][len(campaign['remote_root'])+1:]
  rows[cid]=row;seqs.add(key)
assert rows and len(rows)<=campaign['new_quota']
ordered=sorted(rows.values(),key=lambda r:(-r['ipsae_min'],r['design_id']))
out=args.output/'scored_candidates.jsonl';out.write_bytes(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in ordered).encode())
fields=['design_id','framework_id','structure_method','stratum','parent_id','sequence','length','ipsae_min','monomer_plddt','pose_dockq','computational_quality_pass','final_rank_eligible','source_audit_report']
with (args.output/'scored_candidates.csv').open('w',encoding='utf8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore',lineterminator='\n');w.writeheader();w.writerows(ordered)
passed=[r for r in ordered if r['computational_quality_pass']]
for f in ('progress.json','partition_recovery_20261107.tsv','campaign_finished.json','slurm_cost.json','submission_ids.tsv'):
 if (args.reports/f).exists():(args.output/('source_'+f)).write_bytes((args.reports/f).read_bytes())
scope='partial snapshot of audited v1 scores, not final submission or experimental validation'
finish=args.reports/'campaign_finished.json'
if finish.exists():
 finished=json.loads(finish.read_text())
 if finished['target_reached']:
  progress=finished['progress']
  assert len(rows)==progress['new_reserved']==progress['new_scored_and_audited']==campaign['new_quota']
  assert len(rows)+progress['baseline_scored']==progress['cumulative_scored']==campaign['target']
  assert len(passed)==progress['new_quality_pass'] and not errors
  scope='complete generation and v1 scoring campaign; final candidate selection pending'
write('summary.json',{'scope':scope,'protocol_id':'TNFR1_VHH_ESMFold2_single_seed_v1_20261103','target_range_inclusive':[1,211],'assay_core_mature_rescoring_completed':False,'baseline_scored':1127,'new_scored':len(rows),'cumulative_scored':1127+len(rows),'new_screen_pass':len(passed),'new_passing_parent_ids':sorted({r['parent_id'] for r in passed}),'new_parent_count':len({r['parent_id'] for r in ordered}),'new_scored_by_method':dict(Counter(r['structure_method'] for r in ordered)),'new_scored_by_framework':dict(Counter(r['framework_id'] for r in ordered)),'best_new_ipsae_min':ordered[0]['ipsae_min'],'final_eligible':0,'audit_errors':errors,'source_report_sha256':sources,'coordinate_files_included':False,'output_sha256':{f.name:digest(f) for f in args.output.iterdir() if f.is_file() and f.name not in ('README.md','summary.json')}})
print((args.output/'summary.json').read_text(encoding='utf8').split('  "source_report_sha256"')[0])
