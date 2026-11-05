"""Write explicit, reproducible batch configurations without running any models."""
import json
from pathlib import Path, PurePosixPath

LOCAL = Path(__file__).resolve().parents[1]
ROOT = PurePosixPath('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
folder = LOCAL/'inputs/bulk1'
folder.mkdir(parents=True, exist_ok=True)
rfd, rf, bc = [], [], []
for patch in ('crd2','crd3'):
    for length in ('compact','long'):
        for shard in range(10):
            ident = f'rfd4_{patch}_{length}_{shard:02d}'
            rfd.append(dict(batch_id=ident,method='rfd4',stratum=f'{patch}_{length}',
                            patch=patch,cdr_length=length,seed=110000+len(rfd),
                            n_backbones=20,n_sequences_max=20,
                            input_cif=str(ROOT/'inputs/rfd4'/f'tnfr1_{patch}_{length}.cif'),
                            output_dir=str(ROOT/'generation/bulk1'/ident)))
for patch in ('crd2','crd3'):
    for shard in range(5):
        ident=f'rfantibody_{patch}_{shard:02d}'
        rf.append(dict(batch_id=ident,method='rfantibody',stratum=patch,patch=patch,
                       seed=210000+len(rf)*10,n_backbones=10,n_sequences_max=40,
                       mpnn_seed=42,rf2_seed=310000+len(rf),
                       output_dir=str(ROOT/'generation/bulk1'/ident)))
for patch in ('crd2','crd3'):
    base=json.loads((LOCAL/'inputs'/f'bc2_canary_v2_{patch}.json').read_text())
    for shard in range(20):
        ident=f'bc2_{patch}_{shard:02d}'
        config=json.loads(json.dumps(base))
        config.update(campaign_seed=410000+len(bc),max_trajectories=1,number_of_final_designs=1,
                      project_folder=str(ROOT/'generation/bulk1'/ident))
        (folder/f'{ident}.json').write_text(json.dumps(config,indent=2)+'\n')
        bc.append(dict(batch_id=ident,method='bc2_v2',stratum=patch,patch=patch,
                       seed=config['campaign_seed'],n_backbones=1,n_sequences_max=10,
                       config=str(ROOT/'inputs/bulk1'/f'{ident}.json'),
                       output_dir=config['project_folder']))
# Interleave both epitopes, lengths and methods to spread the first starts.
rfd=sorted(rfd,key=lambda b:(int(b['batch_id'][-2:]),b['stratum']))
rf=sorted(rf,key=lambda b:(int(b['batch_id'][-2:]),b['stratum']))
bc=sorted(bc,key=lambda b:(int(b['batch_id'][-2:]),b['stratum']))
batches=[]
while rfd or rf or bc:
    for group in (rfd,rf,bc,rfd,bc):
        if group: batches.append(group.pop(0))
for task,b in enumerate(batches):
    b['task_id']=task
    b['generation_partition']='a800' if task%3==0 else 'a100'
    b['scoring_partition']='a800' if task%2==0 else 'a100'
campaign=dict(campaign='bulk1',target='TNFR1',root=str(ROOT),
              generation_gpu_concurrency=6,scoring_gpu_concurrency=2,
              generation_walltime='02:00:00',scoring_walltime='01:00:00',
              cpu_walltime='00:30:00',gpu_hours_budget=80,
              expected_attempts=940,sequence_upper_bound=1600,batches=batches)
assert len(batches)==90 and sum(b['n_sequences_max'] for b in batches)==1600
assert sum(b['n_backbones'] for b in batches)==940
assert len({(b['method'],b['seed']) for b in batches})==90
assert all(b['output_dir'].startswith(str(ROOT)+'/generation/bulk1/') and '\\' not in b['output_dir'] for b in batches)
(folder/'campaign.json').write_text(json.dumps(campaign,indent=2)+'\n')
(folder/'batches.tsv').write_bytes(''.join(f"{b['task_id']}\t{b['method']}\t{b['batch_id']}\t{b['generation_partition']}\t{b['scoring_partition']}\n" for b in batches).encode())
for path in (LOCAL/'scripts').iterdir():
    if path.suffix in ('.py','.sh','.slurm'):
        path.write_bytes(path.read_bytes().replace(b'\r\n',b'\n'))
print(json.dumps({'batches':90,'attempts':940,'sequence_upper_bound':1600}))
