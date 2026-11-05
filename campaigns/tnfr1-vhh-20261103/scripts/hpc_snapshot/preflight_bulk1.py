"""Validate frozen campaign inputs using native BC2 settings on a CPU compute node."""
import hashlib
import json
import os
from pathlib import Path
from bindcraft.settings import read_settings, build_design_settings
from bindcraft.preflight import preflight_campaign
from bindcraft.model_weights import model_weights

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
c=json.loads((ROOT/'inputs/bulk1/campaign.json').read_text())
assert len(c['batches'])==90 and c['sequence_upper_bound']==1600
assert sum(b['n_backbones'] for b in c['batches'])==940
assert c['generation_gpu_concurrency']+c['scoring_gpu_concurrency']==8
assert all(b['task_id']==i for i,b in enumerate(c['batches']))
assert len({b['output_dir'] for b in c['batches']})==90
assert len({(b['method'],b['seed']) for b in c['batches']})==90
assert all(Path(b['output_dir']).is_relative_to(ROOT/'generation/bulk1') for b in c['batches'])
p=json.loads((ROOT/'scoring_protocol_v1.json').read_text())
assert p['predictor']['seed']==0 and p['predictor']['num_diffusion_samples']==1
assert p['predictor']['templates'] is False and p['predictor']['target_msa_only']
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
assert sha(ROOT/p['target']['msa'])==p['target']['msa_sha256']
files=[ROOT/'scoring_protocol_v1.json',ROOT/'inputs/bulk1/campaign.json',ROOT/p['target']['msa']]
for input_cif in sorted({b['input_cif'] for b in c['batches'] if b['method']=='rfd4'}):
    assert Path(input_cif).stat().st_size>1000
    files.append(Path(input_cif))
assert len({b['input_cif'] for b in c['batches'] if b['method']=='rfd4'})==4
rfroot=Path('/dssg/home/acct-clshzh/jasoncao1280/RFantibody')
assert '--extra' in (rfroot/'src/rfantibody/cli/inference.py').read_text()
assert 'make_deterministic(i_des)' in (rfroot/'scripts/rfdiffusion_inference.py').read_text()
assert 'design_startnum' in (rfroot/'scripts/config/inference/base.yaml').read_text()
files.extend([rfroot/'src/rfantibody/cli/inference.py',rfroot/'scripts/rfdiffusion_inference.py'])
af2,mpnn=model_weights(download=False)
native_checks=[]
for patch in ('crd2','crd3'):
    group=[b for b in c['batches'] if b['method']=='bc2_v2' and b['patch']==patch]
    assert len(group)==20
    for b in group:
        cfg=json.loads(Path(b['config']).read_text())
        assert cfg['max_trajectories']==cfg['number_of_final_designs']==1
        assert cfg['campaign_seed']==b['seed'] and cfg['project_folder']==b['output_dir']
        assert cfg['min_scaffold_sequence_retained_final']==1.0
        assert cfg['aa_bias']['C']==0 and cfg['losses']['binder_coldspot']['params']['designed_only'] is False
        files.append(Path(b['config']))
    settings=read_settings(group[0]['config'])
    checkpoints=preflight_campaign(settings,settings['project_folder'],af2,mpnn)
    design=build_design_settings(settings)
    assert 'C' in design.binder.omitted_amino_acids
    native_checks.append({'patch':patch,'checkpoints':list(checkpoints)})
# Register previous candidate sequences so bulk scoring does not repeat them.
registry=ROOT/'candidates/registries/bulk1';registry.mkdir(parents=True,exist_ok=True)
for path in (ROOT/'reports').glob('canary_collection_*.json'):
    if path.name.startswith('canary_collection_bulk1_'):continue
    for row in json.loads(path.read_text()).get('candidates',[]):
        digest=hashlib.sha256(row['sequence'].encode()).hexdigest()
        claim=registry/(digest+'.json')
        if not claim.exists():
            with claim.open('x') as f:json.dump({'design_id':row['design_id'],'source_file':row['source_file'],'sequence_sha256':digest,'prior_canary':True},f)
files.extend(sorted((ROOT/'scripts').glob('*bulk1*')))
files.extend([ROOT/'scripts/collect_canary.py',ROOT/'scripts/audit_canary_scores.py',ROOT/'scripts/predict_esmfold2.py'])
report={'ok':True,'job_id':os.environ['SLURM_JOB_ID'],'campaign':'bulk1',
        'native_bc2_checks':native_checks,'hashes':{str(path):sha(path) for path in files if path.is_file()},
        'generator_counts':{name:sum(b['method']==name for b in c['batches']) for name in ('rfd4','rfantibody','bc2_v2')}}
(ROOT/'reports/bulk1_preflight.json').write_text(json.dumps(report,indent=2))
print('BULK1_PREFLIGHT_OK',json.dumps({'batches':90,'sequence_upper_bound':1600}),flush=True)
