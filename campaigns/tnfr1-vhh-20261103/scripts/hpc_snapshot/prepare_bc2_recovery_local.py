import ast
import json
from pathlib import Path,PurePosixPath

LOCAL=Path(__file__).resolve().parents[1]
ROOT=PurePosixPath('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
folder=LOCAL/'inputs/bc2_recovery';folder.mkdir(parents=True,exist_ok=True)
batches=[]
for i in range(40):
    patch=('crd2','crd3')[i%2];shard=i//2
    ident=f'bc2_gpu_{patch}_{shard:02d}'
    config=json.loads((LOCAL/'inputs'/f'bc2_canary_v2_{patch}.json').read_text())
    seed=410000+shard+(20 if patch=='crd3' else 0)
    config.update(campaign_seed=seed,project_folder=str(ROOT/'generation/bc2_recovery'/ident),
                  save_design_trajectory=True,save_design_sequences=True,sparse_output=False)
    (folder/f'{ident}.json').write_text(json.dumps(config,indent=2)+'\n')
    batches.append(dict(task_id=i,batch_id=ident,method='bc2_v2',patch=patch,stratum=patch,seed=seed,
                        output_dir=config['project_folder'],config=str(ROOT/'inputs/bc2_recovery'/f'{ident}.json')))
campaign=dict(campaign='bc2_recovery',batches=batches,baseline_allocated_gpu_hours_ceiling=68.0,
              combined_gpu_hours_limit=80.0,new_gpu_stage_stop_threshold=76.0,
              generation_concurrency=4,scoring_concurrency=2,
              note='Includes saved gradient endpoint before ProteinMPNN plus refolded variants; no change to generation filters or frozen ESMFold2 selection thresholds.')
(folder/'campaign.json').write_text(json.dumps(campaign,indent=2)+'\n')
for path in (LOCAL/'scripts').iterdir():
    if path.name in ('collect_canary.py','audit_canary_scores.py') or 'bc2_recovery' in path.name:
        path.write_bytes(path.read_bytes().replace(b'\r\n',b'\n'))
        if path.suffix=='.py':ast.parse(path.read_text())
assert len({b['output_dir'] for b in batches})==40
print('Prepared40 BC2 GPU-only trajectories with saved gradient endpoints and an80 GPU-hour combined cap.')
