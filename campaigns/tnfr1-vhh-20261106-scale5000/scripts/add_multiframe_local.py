"""Append 18 framework pilots while keeping the cumulative sequence quota."""
import json
from pathlib import Path,PurePosixPath
HERE=Path(__file__).resolve().parents[1];path=HERE/'config/campaign.json';c=json.loads(path.read_text())
assert len(c['batches'])==306,'Already extended or unexpected plan'
ROOT=PurePosixPath(c['remote_root']);batches=[]
for p in ('3K1K','3P0G','1ZVH'):
 route='bc2_native_'+p;first=306+len(batches);c['canary_tasks'][route]=first
 for patch in ('crd2','crd3'):
  for rep in range(3):
   i=306+len(batches);ident=f'bc2_native_{p}_{patch}_{rep:02d}'
   batches.append(dict(task_id=i,batch_id=ident,method='bc2_v2',patch=patch,stratum=patch,framework_id='native_'+p,framework_contract=str(ROOT/'inputs/multiframe'/f'{p}.json'),canary_route=route,canary_task=first,seed=2000000+i,n_sequences_max=11,n_backbones=1,config=str(ROOT/'inputs/multiframe'/f'{ident}.json'),output_dir=str(ROOT/'generation'/ident),generation_partition='a800',scoring_partition='a800'))
c['batches']+=batches;c.update(multiframe_pilot_tasks=[306,323],generation_concurrency_primary=4,generation_concurrency_pilot=2,scoring_concurrency_primary=1,scoring_concurrency_pilot=1)
for p,obj in [(path,c),(HERE/'config/multiframe_batches.json',batches)]:p.write_bytes((json.dumps(obj,indent=2)+'\n').encode())
(HERE/'config/multiframe.tsv').write_bytes(''.join(f"{b['task_id']}\t{b['canary_task']}\n" for b in batches).encode())
print('18 pilots appended; total324; quota remains3873 new / 5000 cumulative; total GPU concurrency6+2')
