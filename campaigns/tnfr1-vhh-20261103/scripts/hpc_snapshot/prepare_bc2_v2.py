"""Prepare and validate stricter CDR-only design objectives on a CPU compute node."""
import json
import os
from pathlib import Path
from bindcraft.settings import read_settings, build_design_settings
from bindcraft.preflight import preflight_campaign
from bindcraft.model_weights import model_weights

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
af2,mpnn = model_weights(download=False)
rows=[]
for index,patch in enumerate(('crd2','crd3')):
    config=json.loads((ROOT/'inputs'/('bc2_canary_'+patch+'.json')).read_text())
    config.update(aa_bias={'C':0,'F':0.7,'W':0.3,'Y':0.9},
                  weights_binder_coldspot=1.0,weights_binder_intra_coldspot=0.0,
                  campaign_seed=303+index,
                  project_folder=str(ROOT/'generation'/('bc2_canary_v2_'+patch)))
    config['losses']={'binder_coldspot':{'params':{'designed_only':False}}}
    path=ROOT/'inputs'/('bc2_canary_v2_'+patch+'.json')
    path.write_text(json.dumps(config,indent=2))
    settings=read_settings(path)
    checkpoints=preflight_campaign(settings,settings['project_folder'],af2,mpnn)
    design=build_design_settings(settings)
    assert 'C' in design.binder.omitted_amino_acids
    assert settings['min_scaffold_sequence_retained_final']==1.0
    assert settings['losses']['binder_coldspot']['params']['designed_only'] is False
    rows.append({'patch':patch,'configuration':str(path),'omitted_new_amino_acids':design.binder.omitted_amino_acids,
                 'fixed_framework_retention':1.0,'whole_framework_target_contact_penalty':1.0,
                 'checkpoints':list(checkpoints)})
(ROOT/'reports/bc2_v2_preflight.json').write_text(json.dumps({'ok':True,'configs':rows},indent=2))
print('BC2_V2_PREFLIGHT_OK',json.dumps(rows),flush=True)
