"""Check the actual installed BC2's campaign interpretation and model files."""
import json
import os
from pathlib import Path
from bindcraft.settings import read_settings, build_design_settings
from bindcraft.preflight import preflight_campaign
from bindcraft.model_weights import model_weights

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
af2, mpnn = model_weights(download=False)
rows = []
for patch in ('crd2','crd3'):
    settings = read_settings(ROOT/'inputs'/('bc2_canary_'+patch+'.json'))
    checkpoints = preflight_campaign(settings, settings['project_folder'], af2, mpnn)
    design = build_design_settings(settings)
    rows.append({'patch': patch, 'checkpoints': list(checkpoints),
                 'scaffold': design.binder.scaffold, 'scaffold_edits': design.binder.scaffold_edits,
                 'binder_chains': list(design.binder_chains)})
result = {'ok': True, 'job_id': os.environ['SLURM_JOB_ID'], 'scope': 'CPU configuration/input preflight, no design', 'configs': rows}
(ROOT/'reports/bc2_input_preflight.json').write_text(json.dumps(result, indent=2))
print('BC2_INPUT_PREFLIGHT_OK',json.dumps(result),flush=True)
