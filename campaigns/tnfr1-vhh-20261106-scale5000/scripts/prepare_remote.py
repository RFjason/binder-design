"""Prepare immutable inputs and backbone-only parent poses on a CPU Slurm node."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import numpy as np
from biotite.structure.io.pdb import PDBFile
from s5k_common import atomic_json,read_json
assert os.environ.get('SLURM_JOB_ID')
ROOT=Path(__file__).resolve().parents[1];c=read_json(ROOT/'config/campaign.json');OLD=Path(c['baseline_root'])
assert len(c['batches'])==306 and c['new_quota']==3873
for sub in ('reports','state','predictions','candidates','cache','generation','logs','inputs/rfd4'):(ROOT/sub).mkdir(parents=True,exist_ok=True)
for name in ('P19438.json','scoring_protocol_v1.json'):shutil.copy2(OLD/name,ROOT/name)
files=['TNFR1_competition_M1_T211.a3m','TNFR1_mature_canonical.pdb','h-NbBCII10_complete.cif','hotspots_crd2.txt','hotspots_crd3.txt']
files += ['rfd4/'+p.name for p in (OLD/'inputs/rfd4').glob('*.cif')]
for name in files:shutil.copy2(OLD/'inputs'/name,ROOT/'inputs'/name)
p=read_json(ROOT/'scoring_protocol_v1.json')
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
assert sha(ROOT/p['target']['msa'])==p['target']['msa_sha256']
registry=read_json(ROOT/'inputs/baseline_registry.json');assert len(registry)==1127
if not (ROOT/'state/registry.json').exists():atomic_json(ROOT/'state/registry.json',registry)
core=['QVQLVESGGGLVQPGGSLRLSC','WFRQAPGQGLEAVA','YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC','WGQGTLVTVS']
parents=read_json(ROOT/'inputs/parents.json');out=[]
for i,parent in enumerate(parents):
 atoms=PDBFile.read(parent['reference_pdb']).get_structure(model=1)
 atoms=atoms[np.isin(atoms.atom_name,['N','CA','C','O'])]
 h=atoms[atoms.chain_id=='H'].copy();t=atoms[atoms.chain_id=='T'].copy()
 assert len(h[h.atom_name=='CA'])==parent['length'] and len(t[t.atom_name=='CA'])==182
 for chain in (h,t):
  old=chain.res_id.copy()
  for j,r in enumerate(chain.res_id[chain.atom_name=='CA']):chain.res_id[old==r]=j+1
 h.set_annotation('b_factor',np.zeros(len(h)));t.set_annotation('b_factor',np.zeros(len(t)))
 pose=h+t;folder=ROOT/'inputs/parents'/f'parent{i:02d}';folder.mkdir(parents=True,exist_ok=True)
 path=folder/'parent.pdb';pdb=PDBFile();pdb.set_structure(pose);pdb.write(path)
 match=re.fullmatch('([ACDEFGHIKLMNPQRSTVWY]{5,30})'.join(map(re.escape,core)),parent['sequence']);assert match
 labels=''.join(f'REMARK PDBinfo-LABEL: {j+1:4d} H{g}\n' for g in (1,2,3) for j in range(match.start(g),match.end(g)))
 with path.open('a') as f:f.write(labels)
 out.append(dict(index=i,design_id=parent['design_id'],source_sha256=sha(parent['reference_pdb']),input_sha256=sha(path),loops=parent['cdr_lengths']))
atomic_json(ROOT/'reports/prepared.json',{'ok':True,'job_id':os.environ['SLURM_JOB_ID'],'parents':out,'msa_sha256':sha(ROOT/p['target']['msa'])})
print('S5K_INPUTS_PREPARED',flush=True)
