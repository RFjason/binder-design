"""Independent IMGT numbering of the 19 originals and their tail-repaired variants."""
import csv
import json
import os
from pathlib import Path
import sys
from zipfile import ZipFile
assert os.environ.get('SLURM_JOB_ID'),'Slurm CPU compute node only'
ROOT=Path(__file__).resolve().parent
deps=ROOT/'runtime';deps.mkdir(exist_ok=True)
for path in (ROOT/'wheels').glob('*.whl'):
 with ZipFile(path) as z:
  for name in z.namelist():
   assert (deps/name).resolve().is_relative_to(deps.resolve())
  z.extractall(deps)
sys.path.insert(0,str(deps))
from anarcii import Anarcii
items=list(csv.DictReader((ROOT/'candidates.csv').open()))
seqs=[(r['design_id'],r['sequence']) for r in items]
seqs += [(r['design_id']+'_fr4S',r['sequence']+'S') for r in items]
model=Anarcii(seq_type='antibody',cpu=True,ncpu=2,batch_size=32,verbose=False)
result=model.number(seqs)
(ROOT/'anarcii_raw.json').write_text(json.dumps(result,default=str,indent=2))
rows=[]
for sid,seq in seqs:
 rec=result[sid];assert rec.get('numbering'),(sid,rec)
 numbered=rec['numbering'];real=[(p,aa) for p,aa in numbered if aa!='-']
 cdrs={name:''.join(aa for (p,ins),aa in numbered if lo<=p<=hi and aa!='-') for name,lo,hi in [('CDR1',27,38),('CDR2',56,65),('CDR3',105,117)]}
 rows.append({'design_id':sid,'sequence':seq,'last_imgt':real[-1][0],'last_aa':real[-1][1],'cdrs':cdrs,'score':rec.get('score'),'error':rec.get('error')})
original=[r for r in rows if not r['design_id'].endswith('_fr4S')]
groups={}
for r in original:
 key=tuple(r['cdrs'][k] for k in ('CDR1','CDR2','CDR3'));groups.setdefault(key,[]).append(r['design_id'])
result={'job_id':os.environ['SLURM_JOB_ID'],'anarcii_version':'2.0.8','original_count':len(original),'unique_full_CDR_triplets':len(groups),'duplicates':[v for v in groups.values() if len(v)>1],'rows':rows}
(ROOT/'independent_numbering.json').write_text(json.dumps(result,indent=2))
print('INDEPENDENT_IMGT_RESULT',json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
for r in original:print(r['design_id'],r['last_imgt'],r['cdrs'],flush=True)
