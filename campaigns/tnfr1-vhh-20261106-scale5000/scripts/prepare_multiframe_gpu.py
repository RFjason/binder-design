"""Complete the mature source scaffold by an explicit GPU monomer fold."""
import hashlib,json,os,subprocess,sys
from pathlib import Path
import numpy as np
from biotite.structure.io.pdbx import CIFFile,get_structure,set_structure
from s5k_common import atomic_json,read_json
assert os.environ.get('SLURM_JOB_ID') and os.environ.get('CUDA_VISIBLE_DEVICES')
ROOT=Path(__file__).resolve().parents[1]
b=json.loads(sys.argv[1]);c=read_json(b['framework_contract']);cfg=read_json(b['config']);dest=Path(cfg['binder_scaffold'])
if dest.exists():raise SystemExit(0)
ident='framework_'+c['framework_id']
manifest=ROOT/'inputs/multiframe'/(ident+'_fold.json')
atomic_json(manifest,[{'id':ident,'candidate_id':ident,'chains':[{'id':'B','sequence':c['sequence']}]}])
subprocess.run(['bash',str(ROOT/'scripts/predict.slurm'),str(manifest)],check=True,cwd=ROOT)
folder=ROOT/'predictions'/ident;summary=read_json(folder/'summary.json')
assert summary['ok'] and summary['chain_plddt']['B']>=70
atoms=get_structure(CIFFile.read(folder/'model.cif'),model=1)
assert len(atoms[atoms.atom_name=='CA'])==len(c['sequence'])
sg=atoms[(atoms.res_name=='CYS')&(atoms.atom_name=='SG')]
assert len(sg)==2 and np.linalg.norm(sg.coord[0]-sg.coord[1])<=3
atoms.chain_id[:]='H';ca=atoms[atoms.atom_name=='CA'];old=atoms.res_id.copy()
for i,resid in enumerate(ca.res_id):atoms.res_id[old==resid]=i+1
dest.parent.mkdir(parents=True,exist_ok=True);cif=CIFFile();set_structure(cif,atoms);cif.write(dest)
atomic_json(dest.with_suffix('.provenance.json'),{'framework_id':c['framework_id'],'source_contract':b['framework_contract'],'source_sequence_sha256':hashlib.sha256(c['sequence'].encode()).hexdigest(),'model_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'fold_summary':summary,'coordinate_source':'ESMFold2 monomer; full mature sequence; not original PDB coordinates'})
