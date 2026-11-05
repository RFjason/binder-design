"""Prepare exact-sequence prediction cases and normalize text upload files."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
target = json.loads((ROOT / 'sources/P19438.json').read_text())['sequence']['value'][:211]
msa_path = ROOT / 'inputs/TNFR1_competition_M1_T211.a3m'
lines = msa_path.read_text().splitlines()
trimmed = []
for line in lines:
    if line.startswith('>'):
        trimmed.append(line)
        continue
    position, kept = 0, []
    for letter in line:
        if letter.isupper() or letter == '-':
            position += 1
        if position > 29:
            kept.append(letter)
    assert position == 211
    trimmed.append(''.join(kept))
assert trimmed[1] == target[29:]
(ROOT / 'inputs/TNFR1_mature_L30_T211.a3m').write_text('\n'.join(trimmed) + '\n')
(ROOT / 'inputs/TNFR1_mature_L30_T211.fasta').write_text('>TNFR1_mature_L30_T211\n' + target[29:] + '\n')
cases = [
    {'id': 'target_mature', 'purpose': 'full mature ectodomain design context',
     'chains': [{'id': 'A', 'sequence': target[29:], 'msa': 'inputs/TNFR1_mature_L30_T211.a3m'}]},
    {'id': 'target_competition', 'purpose': 'competition input diagnostics',
     'chains': [{'id': 'A', 'sequence': target, 'msa': 'inputs/TNFR1_competition_M1_T211.a3m'}]},
]
(ROOT / 'inputs/target_prediction_cases.json').write_text(json.dumps(cases, indent=2))
for path in (ROOT / 'scripts').iterdir():
    if path.suffix in ('.py', '.slurm', '.sh'):
        path.write_bytes(path.read_bytes().replace(b'\r\n', b'\n'))
manifest = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for folder in ('inputs', 'scripts') for p in (ROOT / folder).iterdir() if p.is_file()}
(ROOT / 'reports/initial_upload_sha256.json').write_text(json.dumps(manifest, indent=2))
print(json.dumps({'target_length': len(target), 'mature_length': len(target[29:]), 'files': len(manifest)}))
