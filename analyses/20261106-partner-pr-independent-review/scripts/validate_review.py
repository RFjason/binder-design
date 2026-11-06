"""Verify the evidence package and recompute its claims using the standard library."""
import csv
from difflib import SequenceMatcher
import hashlib
from itertools import combinations
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
def read(path):return json.loads(path.read_text(encoding='utf8'))
def table(name):return list(csv.DictReader((ROOT/'results'/name).open(encoding='utf8',newline='')))
manifest=read(ROOT/'manifest.json')
for name,digest in manifest['files'].items():
 path=ROOT/name
 assert path.resolve().is_relative_to(ROOT.resolve())
 assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,name
for line in (ROOT/'sources/pr1/analyses/20261105-compliance-novelty-audit/input_hashes.txt').read_text().splitlines():
 digest,name=line.split(None,1)
 path=ROOT/'sources/pr1'/name if name.startswith('analyses/') else REPO/name
 assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,name

base=[json.loads(s) for s in (REPO/'campaigns/tnfr1-vhh-20261103/results/all_candidates.jsonl').read_text().splitlines()]
short={r['design_id']:r for r in base if r['computational_quality_pass']}
assert len(base)==1127 and len(short)==19
numbered=read(ROOT/'results/independent_numbering.json');raw=read(ROOT/'results/anarcii_raw.json')
rows=numbered['rows'];assert len(rows)==38
original=[r for r in rows if not r['design_id'].endswith('_fr4S')]
assert len(original)==19 and len({tuple(r['cdrs'][k] for k in ('CDR1','CDR2','CDR3')) for r in original})==19
for r in rows:
 repaired=r['design_id'].endswith('_fr4S');sid=r['design_id'][:-5] if repaired else r['design_id']
 assert r['sequence']==short[sid]['sequence']+('S' if repaired else '')
 assert r['last_imgt'][0]==(128 if repaired else 127)
 numbering=raw[r['design_id']]['numbering']
 for key,lo,hi in [('CDR1',27,38),('CDR2',56,65),('CDR3',105,117)]:
  seq=''.join(aa for (pos,ins),aa in numbering if lo<=pos<=hi and aa!='-')
  assert r['cdrs'][key]==seq
for x,y in combinations(original,2):
 a=''.join(x['cdrs'][k] for k in ('CDR1','CDR2','CDR3'));b=''.join(y['cdrs'][k] for k in ('CDR1','CDR2','CDR3'))
 match=SequenceMatcher(None,a,b,autojunk=False).find_longest_match(0,len(a),0,len(b))
 assert match.size/max(len(a),len(b))<=.9
for x,y in [('de38b5b8dffd','203250d5c341'),('9a7253825d68','bf4208b6c1a3'),('60371bf2a410','14c12a1d7ac8')]:
 a=next(r for r in original if r['design_id'].endswith(x));b=next(r for r in original if r['design_id'].endswith(y))
 assert a['cdrs']['CDR1']==b['cdrs']['CDR1'] and a['cdrs']['CDR3']==b['cdrs']['CDR3']
 assert a['cdrs']['CDR2']!=b['cdrs']['CDR2']

cal=table('calibration_seed0_per_design.csv');assert len(cal)==1320 and len({r['uuid'] for r in cal})==1320
assert all(r['seed']=='0' for r in cal)
assert sum(r['binder_final']=='True' for r in cal)==354
for row in table('threshold_calibration.csv'):
 threshold=float(row['threshold']);counts=dict(TP=0,FP=0,FN=0,TN=0)
 for r in cal:
  positive=r['binder_final']=='True';passes=float(r['ipsae_min'])>=threshold
  key=('TP' if positive else 'FP') if passes else ('FN' if positive else 'TN');counts[key]+=1
 assert all(counts[k]==int(row[k]) for k in counts)
changed=[r for r in cal if float(r['ipsae_min'])<.6<=float(r['five_seed_max'])]
assert len(changed)==73 and sum(r['binder_final']=='True' for r in changed)==13

for method in ('rfd4','rfantibody','bindcraft2'):
 group=[r for r in base if r['structure_method']==method]
 assert sum(r['ipsae_min']>=.6 for r in group)==sum(r['computational_quality_pass'] for r in group)
for patch in ('crd2','crd3'):
 group=[r for r in base if r['structure_method']=='bindcraft2' and r['stratum']==patch and r['computational_quality_pass']]
 assert len({r['parent_id'] for r in group})==2
assert read(ROOT/'results/local_checks.json')['missing_novelty_backend_exit_code']==0
assert table('novelty_backend_repro.csv')[0]['结果']=='PASS'
print(f'PASS: {len(manifest["files"])} evidence hashes; 38 raw numberings; 19 distinct CDR triplets; 1320 calibration rows; threshold tables; route/parent claims.')
