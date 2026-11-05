"""Local source retrieval for planning only; never submits Slurm jobs."""
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT.parent
VENDOR = WORK / 'IL7RA_20261102' / '.vendor'
sys.path.insert(0, str(VENDOR))
import requests

SKILLS = Path('C:/Users/caoji/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills')
for directory in ('sources', 'reports', 'inputs'):
    (ROOT / directory).mkdir(parents=True, exist_ok=True)

def api(skill, request, filename):
    request.update(save_raw=True, raw_output_path=str(ROOT / 'sources' / filename), max_depth=1)
    result = subprocess.run(
        [sys.executable, str(SKILLS / skill / 'scripts/rest_request.py')],
        input=json.dumps(request), capture_output=True, text=True,
        env={**os.environ, 'PYTHONPATH': str(VENDOR)}, timeout=70)
    (ROOT / 'reports' / (filename + '.request.json')).write_text(result.stdout, encoding='utf-8')
    response = json.loads(result.stdout)
    if result.returncode or not response.get('ok'):
        raise RuntimeError(result.stdout + result.stderr)

api('uniprot-skill', {'base_url': 'https://rest.uniprot.org', 'path': 'uniprotkb/P19438',
                     'params': {'format': 'json'}}, 'P19438.json')
api('rcsb-pdb-skill', {'base_url': 'https://data.rcsb.org/rest/v1', 'path': 'core/entry/9S3I'}, '9S3I_entry.json')
api('rcsb-pdb-skill', {'base_url': 'https://data.rcsb.org/rest/v1', 'path': 'core/polymer_entity/9S3I/1'}, '9S3I_entity1.json')
downloads = {
    '9S3I.pdb': 'https://files.rcsb.org/download/9S3I.pdb',
    '1TNR.pdb': 'https://files.rcsb.org/download/1TNR.pdb',
    'Anthropic_single_target_EGFR_prompt.md': 'https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/resolve/main/prompts/prompts/single_target/EGFR.md',
}
for filename, url in downloads.items():
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    (ROOT / 'sources' / filename).write_bytes(response.content)

record = json.loads((ROOT / 'sources/P19438.json').read_text(encoding='utf-8'))
sequence = record['sequence']['value']
target = sequence[:211]
assert len(target) == 211 and target[0] == 'M' and target[-1] == 'T'
(ROOT / 'inputs/TNFR1_competition_M1_T211.fasta').write_text('>TNFR1_P19438_1_M1_T211\n' + target + '\n', encoding='ascii')
summary = {'uniprot': record['primaryAccession'], 'full_length': len(sequence),
           'competition_length': len(target), 'features': [f for f in record.get('features', [])
           if f['type'] in ('Signal', 'Topological domain', 'Domain', 'Glycosylation', 'Disulfide bond')],
           'status': 'planning_only_no_compute_submitted'}
(ROOT / 'reports/target_planning_audit.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary, indent=2))
