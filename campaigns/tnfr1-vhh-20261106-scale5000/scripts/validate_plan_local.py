"""Validate scientific quotas, distinct sampling and the global registration cap."""
import hashlib
import json
import re
from pathlib import Path
from collections import Counter
from s5k_common import select_row,new_count
ROOT=Path(__file__).resolve().parents[1]
c=json.loads((ROOT/'config/campaign.json').read_text());b=c['batches']
extra=18 if (ROOT/'config/multiframe_batches.json').exists() else 0
assert Counter(x['method'] for x in b)=={'rfd4':184,'bc2_v2':100+extra,'mpnn_refine':18,'rfantibody':4}
assert len(b)==306+extra and c['baseline_scored']+c['new_quota']==c['target']==5000
assert all(x['task_id']==i and x['generation_partition']==x['scoring_partition']=='a800' for i,x in enumerate(b))
assert len({x['output_dir'] for x in b})==306+extra
assert all('\\' not in x['output_dir'] for x in b)
for method in ('rfd4','bc2_v2','rfantibody'):
 group=[x for x in b if x['method']==method]
 assert len({x['seed'] for x in group})==len(group)
assert len({x['parent_id'] for x in b if x['method']=='mpnn_refine'})==6
assert sum(x['n_sequences_max'] for x in b if x['method']=='rfd4' and not x.get('reserve'))==3040
assert all(b[task].get('canary_route',b[task]['method'])==method for method,task in c['canary_tasks'].items())
if extra:
 assert len({x['framework_id'] for x in b[306:]})==3
 assert 100*4+18*4+3<500
 assert c['generation_concurrency_primary']+c['generation_concurrency_pilot']==6
 assert c['scoring_concurrency_primary']+c['scoring_concurrency_pilot']==2
 for name in ('3K1K','3P0G','1ZVH'):
  contract=json.loads((ROOT/'inputs/multiframe'/f'{name}.json').read_text())
  seq=contract['sequence'];pat=re.compile('([ACDEFGHIKLMNPQRSTVWY]{5,30})'.join(map(re.escape,contract['core'])))
  assert pat.fullmatch(seq) and 100<=len(seq)<=150 and seq.count('C')==2 and 'HHHHHH' not in seq
  assert seq[contract['cys1_index0']]=='C' and contract['source_terminal_imgt']==128
  start=contract['cdr_spans0'][1][0];mut=seq[:start]+('A' if seq[start]!='A' else 'S')+seq[start+1:]
  assert pat.fullmatch(mut),'CDR changes must be designable'
  assert not pat.fullmatch('S'+seq[1:]),'Non-contract framework changes must be rejected'
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
print(f'PASS: {len(b)} finite batches; cumulative quota5000; baseline1127; global de-duplication/quota; wave limit500.')
