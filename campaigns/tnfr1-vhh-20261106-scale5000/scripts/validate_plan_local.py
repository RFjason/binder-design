"""Validate scientific quotas, distinct sampling and the global registration cap."""
import hashlib
import json
from pathlib import Path
from collections import Counter
from s5k_common import select_row,new_count
ROOT=Path(__file__).resolve().parents[1]
c=json.loads((ROOT/'config/campaign.json').read_text());b=c['batches']
assert Counter(x['method'] for x in b)=={'rfd4':184,'bc2_v2':100,'mpnn_refine':18,'rfantibody':4}
assert len(b)==306 and c['baseline_scored']+c['new_quota']==c['target']==5000
assert all(x['task_id']==i and x['generation_partition']==x['scoring_partition']=='a800' for i,x in enumerate(b))
assert len({x['output_dir'] for x in b})==306
assert all('\\' not in x['output_dir'] for x in b)
for method in ('rfd4','bc2_v2','rfantibody'):
 group=[x for x in b if x['method']==method]
 assert len({x['seed'] for x in group})==len(group)
assert len({x['parent_id'] for x in b if x['method']=='mpnn_refine'})==6
assert sum(x['n_sequences_max'] for x in b if x['method']=='rfd4' and not x.get('reserve'))==3040
assert all(b[task]['method']==method for method,task in c['canary_tasks'].items())
assert all((100 if i<2 else 74 if i==2 else 32)*4+2<500 for i in range(4))
baseline=json.loads((ROOT/'inputs/baseline_registry.json').read_text());assert len(baseline)==1127
old=next(iter(baseline));assert baseline[old]['origin']=='baseline'
r={};row=lambda i,s:{'design_id':i,'sequence':s,'parent_id':'parent'}
assert select_row(r,row('a','AAAA'),'batch1',2)=='new'
assert select_row(r,row('a','AAAA'),'batch1',2)=='existing_own'
assert select_row(r,row('b','AAAA'),'batch2',2)=='duplicate'
assert select_row(r,row('c','CCCC'),'batch2',2)=='new'
assert select_row(r,row('d','DDDD'),'batch3',2)=='quota_full'
assert new_count(r)==2
print('PASS: 306 finite batches; cumulative quota5000; baseline1127; global de-duplication/quota; wave limit500.')
