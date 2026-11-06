"""Create an explicit finite sampling graph locally; no remote computation."""
import hashlib
import json
from pathlib import Path, PurePosixPath

HERE=Path(__file__).resolve().parents[1]
REPO=HERE.parents[1]
OLD_LOCAL=REPO.parent/'TNFR1_VHH_20261103'
OLD='/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103'
ROOT=PurePosixPath('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261106_scale5000')
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes((json.dumps(value,indent=2,ensure_ascii=False)+'\n').encode())
baseline=[json.loads(line) for line in (REPO/'campaigns/tnfr1-vhh-20261103/results/all_candidates.jsonl').read_text().splitlines()]
assert len(baseline)==1127
registry={}
for r in baseline:
    # Normalizing the old truncated tail also prevents counting pure tail repairs.
    seq=r['sequence'];seq+= 'S' if seq.endswith('WGQGTLVTVS') else ''
    key=hashlib.sha256(seq.encode()).hexdigest()
    assert key not in registry
    registry[key]={'origin':'baseline','design_id':r['design_id']}
write(HERE/'inputs/baseline_registry.json',registry)
ids=['TNFR1_bulk1_rfd4_crd2_long_03_04ca5213f78c',
     'TNFR1_bulk1_rfd4_crd2_compact_00_8fced7a5fd5b',
     'TNFR1_bc2_recovery_bc2_gpu_crd3_12_083b9fc2dcc0',
     'TNFR1_bc2_recovery_bc2_gpu_crd2_08_60371bf2a410',
     'TNFR1_bc2_recovery_bc2_gpu_crd2_01_7d878eb23d8e',
     'TNFR1_bulk1_rfd4_crd3_long_01_ab9de9edd9d5']
parents=[next(r for r in baseline if r['design_id']==i) for i in ids]
write(HERE/'inputs/parents.json',parents)
rfd=[];bc=[];mp=[];rf=[];reserve=[]
for shard in range(38):
 for patch in ('crd2','crd3'):
  for length in ('compact','long'):
   ident=f'rfd4_{patch}_{length}_{shard:02d}'
   rfd.append(dict(batch_id=ident,method='rfd4',patch=patch,stratum=f'{patch}_{length}',seed=510000+len(rfd),n_sequences_max=20,n_backbones=20,input_cif=str(ROOT/'inputs/rfd4'/f'tnfr1_{patch}_{length}.cif')))
for shard in range(50):
 for patch in ('crd2','crd3'):
  ident=f'bc2_{patch}_{shard:02d}'
  cfg=json.loads((OLD_LOCAL/'inputs/bc2_recovery'/f'bc2_gpu_{patch}_00.json').read_text().replace(OLD,str(ROOT)))
  cfg.update(campaign_seed=610000+len(bc),project_folder=str(ROOT/'generation'/ident))
  write(HERE/'inputs/bc2'/f'{ident}.json',cfg)
  bc.append(dict(batch_id=ident,method='bc2_v2',patch=patch,stratum=patch,seed=cfg['campaign_seed'],n_sequences_max=11,n_backbones=1,config=str(ROOT/'inputs/bc2'/f'{ident}.json')))
for i,p in enumerate(parents):
 for t in (.1,.2,.4):
  ident=f'mpnn_parent{i:02d}_t{int(t*10):02d}'
  mp.append(dict(batch_id=ident,method='mpnn_refine',stratum=f'parent{i:02d}',seed=42,temperature=t,n_sequences_max=20,n_backbones=0,parent_id=p['parent_id'],parent_design_id=p['design_id'],parent_structure_method=p['structure_method'],parent_input=str(ROOT/'inputs/parents'/f'parent{i:02d}'),parent_index=i))
for i in range(4):
 patch='crd2' if i%2==0 else 'crd3';ident=f'rfantibody_{patch}_{i:02d}'
 rf.append(dict(batch_id=ident,method='rfantibody',patch=patch,stratum=patch,seed=710000+10*i,rf2_seed=810000+i,n_sequences_max=16,n_backbones=4))
for i in range(32):
 patch=('crd2','crd3')[i%2];length=('compact','long')[(i//2)%2];ident=f'reserve_{patch}_{length}_{i:02d}'
 reserve.append(dict(batch_id=ident,method='rfd4',patch=patch,stratum=f'{patch}_{length}',seed=910000+i,n_sequences_max=20,n_backbones=20,input_cif=str(ROOT/'inputs/rfd4'/f'tnfr1_{patch}_{length}.cif'),reserve=True))
# First task of each route is its same-spec canary.
batches=[rfd.pop(0),bc.pop(0),mp.pop(0),rf.pop(0)]
while rfd or bc or mp or rf:
 for group in (rfd,bc,rfd,mp,rf):
  if group:batches.append(group.pop(0))
assert len(batches)==274
batches+=reserve
canaries={b['method']:i for i,b in enumerate(batches[:4])}
for i,b in enumerate(batches):
 b.update(task_id=i,canary_task=canaries[b['method']],output_dir=str(ROOT/'generation'/b['batch_id']),generation_partition='a800',scoring_partition='a800')
assert len(batches)==306
if (HERE/'config/multiframe_batches.json').exists():
 extra=json.loads((HERE/'config/multiframe_batches.json').read_text())
 batches+=extra
 for b in extra:canaries[b['canary_route']]=b['canary_task']
c=dict(campaign='scale5000',remote_root=str(ROOT),baseline_root=OLD,baseline_scored=1127,new_quota=3873,target=5000,batches=batches,primary_batches=274,canary_tasks=canaries,generation_concurrency=6,scoring_concurrency=2,baseline_allocated_gpu_hours=74.19666666666662,total_gpu_hours_budget=200,new_gpu_stage_stop_threshold=193)
if len(batches)>306:c.update(multiframe_pilot_tasks=[306,323],generation_concurrency_primary=4,generation_concurrency_pilot=2,scoring_concurrency_primary=1,scoring_concurrency_pilot=1)
write(HERE/'config/campaign.json',c)
(HERE/'config/batches.tsv').write_bytes(''.join(f"{b['task_id']}\t{b['method']}\t{b['batch_id']}\t{b['generation_partition']}\t{b['scoring_partition']}\t{b['canary_task']}\n" for b in batches).encode())
for name in ('predict_esmfold2.py','predict.slurm','audit_canary_scores.py','run_bc2_gpu_required.py'):
 if (HERE/'scripts'/name).exists():continue
 s=(OLD_LOCAL/'scripts'/name).read_text().replace(OLD,str(ROOT))
 if name=='audit_canary_scores.py':
  s=s.replace("choices=['bulk1','bc2_recovery']","choices=['bulk1','bc2_recovery','scale5000']")
  s=s.replace("'WGQGTLVTVS']","'WGQGTLVTVSS']")
 if name=='predict_esmfold2.py':
  s=s.replace("ROOT = Path(","assert torch.cuda.is_available(), 'CUDA required; no CPU fallback'\nROOT = Path(",1)
 (HERE/'scripts'/name).write_bytes(s.replace('\r\n','\n').encode())
for path in (HERE/'scripts').glob('*'):
 if path.suffix in ('.py','.sh','.slurm'):path.write_bytes(path.read_bytes().replace(b'\r\n',b'\n'))
print(json.dumps({'batches':len(batches),'new_quota':3873,'canary_tasks':canaries,'baseline_unique':len(registry)}))
