"""Create a small diagnostic score and stricter next BC2 canary configurations."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
collection = json.loads((ROOT/'reports/canary_collection_bc2.json').read_text())
target = json.loads((ROOT/'sources/P19438.json').read_text())['sequence']['value'][:211]
candidate = collection['candidates'][0]
cases = [
    {'id': candidate['design_id']+'_complex', 'candidate_id': candidate['design_id'],
     'purpose': 'one rejected candidate, scoring-adapter diagnostic only',
     'chains': [{'id':'A','sequence':target,'msa':'inputs/TNFR1_competition_M1_T211.a3m'},
                {'id':'B','sequence':candidate['sequence']}]},
    {'id':candidate['design_id']+'_monomer', 'candidate_id':candidate['design_id'],
     'purpose':'diagnostic only', 'chains':[{'id':'B','sequence':candidate['sequence']}]}
]
(ROOT/'inputs/canary_scoring_bc2.json').write_text(json.dumps(cases,indent=2))
summary = {'n_generated_unique':len(collection['candidates']),
           'lengths': sorted({row['length'] for row in collection['candidates']}),
           'cysteine_counts': sorted({row['cysteine_count'] for row in collection['candidates']}),
           'cysteine_geometry_pass':sum(row['cysteine_geometry_ok'] for row in collection['candidates']),
           'selected_scoring_diagnostic':candidate['design_id'], 'not_eligible_for_submission':True}
(ROOT/'reports/bc2_v1_progress_summary.json').write_text(json.dumps(summary,indent=2))
for path in (ROOT/'scripts').iterdir():
    if path.suffix in ('.slurm','.sh','.py'):
        path.write_bytes(path.read_bytes().replace(b'\r\n',b'\n'))
print(json.dumps(summary))
