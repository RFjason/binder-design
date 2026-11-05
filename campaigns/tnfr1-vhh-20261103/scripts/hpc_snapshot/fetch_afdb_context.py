"""Retrieve official full-length context; inference and structure analysis stay in Slurm."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
vendor = ROOT.parent / 'IL7RA_20261102/.vendor'
sys.path.insert(0, str(vendor))
import requests

client = Path('C:/Users/caoji/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills/alphafold-skill/scripts/rest_request.py')
request = {'base_url': 'https://alphafold.ebi.ac.uk/api', 'path': 'prediction/P19438',
           'save_raw': True, 'raw_output_path': str(ROOT / 'sources/AFDB_P19438_metadata.json'), 'max_depth': 1}
call = subprocess.run([sys.executable, str(client)], input=json.dumps(request), capture_output=True,
                      text=True, env={**os.environ, 'PYTHONPATH': str(vendor)}, timeout=60)
(ROOT / 'reports/AFDB_request.json').write_text(call.stdout)
assert call.returncode == 0 and json.loads(call.stdout)['ok'], call.stdout
records = json.loads((ROOT / 'sources/AFDB_P19438_metadata.json').read_text())
canonical = json.loads((ROOT / 'sources/P19438.json').read_text())['sequence']['value']
record = next(x for x in records if x.get('uniprotStart') == 1 and x.get('uniprotEnd') == 455)
assert record['sequence'] == canonical
response = requests.get(record['pdbUrl'], timeout=50)
response.raise_for_status()
path = ROOT / 'sources/AFDB_P19438_full.pdb'
path.write_bytes(response.content)
lines = [x for x in response.text.splitlines() if x.startswith('ATOM  ') and 30 <= int(x[22:26]) <= 211]
# Source retains canonical numbering, mature design-input parser handles this explicitly.
(ROOT / 'inputs/TNFR1_AFDB_mature_L30_T211.pdb').write_text('\n'.join(lines) + '\nTER\nEND\n')
print(json.dumps({'source': record['pdbUrl'], 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                  'mature_atom_lines': len(lines), 'usage': 'design context only, not final scoring template'}))
