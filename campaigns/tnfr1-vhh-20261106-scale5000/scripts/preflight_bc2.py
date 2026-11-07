"""Use installed BC2's native configuration validator on a CPU node."""
import hashlib
import json
import os
from pathlib import Path
from bindcraft.settings import read_settings,build_design_settings
from bindcraft.preflight import preflight_campaign
from bindcraft.model_weights import model_weights
assert os.environ.get('SLURM_JOB_ID')
ROOT=Path(__file__).resolve().parents[1]
c=json.loads((ROOT/'config/campaign.json').read_text());af2,mpnn=model_weights(download=False)
rows=[]
for patch in ('crd2','crd3'):
 group=[b for b in c['batches'] if b['method']=='bc2_v2' and b['patch']==patch and not b.get('framework_contract')];assert len(group)==50
 for b in group:
  cfg=json.loads(Path(b['config']).read_text());assert cfg['campaign_seed']==b['seed']
  assert cfg['project_folder']==b['output_dir'] and cfg['max_trajectories']==1
  assert cfg['save_design_trajectory'] and cfg['min_scaffold_sequence_retained_final']==1
 s=read_settings(group[0]['config']);checkpoints=preflight_campaign(s,s['project_folder'],af2,mpnn);d=build_design_settings(s)
 assert 'C' in d.binder.omitted_amino_acids
 rows.append({'patch':patch,'checkpoints':list(checkpoints)})
hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ('config','scripts') for p in (ROOT/folder).glob('*') if p.is_file()}
(ROOT/'reports/preflight.json').write_text(json.dumps({'ok':True,'job_id':os.environ['SLURM_JOB_ID'],'native_bc2':rows,'hashes':hashes},indent=2))
print('S5K_NATIVE_PREFLIGHT_OK',flush=True)
