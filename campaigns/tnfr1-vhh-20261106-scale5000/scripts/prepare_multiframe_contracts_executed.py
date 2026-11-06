"""Resolve complete mature framework contracts using IMGT on a Slurm CPU node."""
import hashlib,json,os,re,sys
from pathlib import Path
assert os.environ.get('SLURM_JOB_ID')
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,'/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_pr_review_20261106/runtime')
from anarcii import Anarcii
records=[]
for pdb in ('3K1K','3P0G','1ZVH'):
 text=(ROOT/'sources'/f'{pdb}.fasta').read_text()
 entries=[x.splitlines() for x in text.split('>')[1:]]
 entry=next(x for x in entries if x[0].startswith(pdb+'_2|'))
 records.append((pdb,''.join(entry[1:])))
res=Anarcii(seq_type='antibody',cpu=True,ncpu=2,batch_size=8,verbose=False).number(records)
contracts=[]
for pdb,raw in records:
 numbered=[(p,aa) for p,aa in res[pdb]['numbering'] if aa!='-' and 1<=p[0]<=128]
 seq=''.join(aa for p,aa in numbered)
 assert numbered[0][0][0]==1 and numbered[-1][0][0]==128,(pdb,numbered[0],numbered[-1])
 assert 100<=len(seq)<=150 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY') and seq.count('C')==2
 indices={p[0]:i for i,(p,aa) in enumerate(numbered) if not p[1].strip()}
 assert seq[indices[23]]==seq[indices[104]]=='C'
 spans=[]
 for lo,hi in ((27,38),(56,65),(105,117)):
  positions=[i for i,(p,aa) in enumerate(numbered) if lo<=p[0]<=hi]
  assert positions and positions==list(range(min(positions),max(positions)+1))
  spans.append([min(positions),max(positions)+1])
 core=[seq[:spans[0][0]],seq[spans[0][1]:spans[1][0]],seq[spans[1][1]:spans[2][0]],seq[spans[2][1]:]]
 c={'framework_id':'native_'+pdb,'pdb':pdb,'source_url':'https://www.rcsb.org/structure/'+pdb,'source_fasta_sha256':hashlib.sha256((ROOT/'sources'/f'{pdb}.fasta').read_bytes()).hexdigest(),'source_pdb_sha256':hashlib.sha256((ROOT/'sources'/f'{pdb}.pdb').read_bytes()).hexdigest(),'raw_entity_sequence':raw,'sequence':seq,'core':core,'cdr_spans0':spans,'cys1_index0':indices[23],'source_terminal_imgt':128,'mature_tail_complete':True,'coordinate_preparation':'seed0 ESMFold2 monomer within the first GPU generation task; not original experimental coordinates','framework_family_claim':'native PDB-derived scaffold; not the engineered sybody library consensus'}
 (ROOT/'contracts').mkdir(exist_ok=True);(ROOT/'contracts'/f'{pdb}.json').write_text(json.dumps(c,indent=2));contracts.append(c)
assert len({tuple(c['core']) for c in contracts})==3
(ROOT/'preflight.json').write_text(json.dumps({'ok':True,'job_id':os.environ['SLURM_JOB_ID'],'contracts':contracts},indent=2))
print('MULTIFRAME_CONTRACTS_OK',[(c['framework_id'],len(c['sequence']),[b-a for a,b in c['cdr_spans0']]) for c in contracts],flush=True)
