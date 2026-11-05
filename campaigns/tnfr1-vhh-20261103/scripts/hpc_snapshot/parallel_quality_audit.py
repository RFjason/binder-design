"""Independent geometry, paratope, liability and surface descriptors on a CPU node."""
import hashlib
import itertools
import json
import os
from pathlib import Path
import re

import numpy as np
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile,get_structure

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
FRAMEWORK=['QVQLVESGGGLVQPGGSLRLSC','WFRQAPGQGLEAVA','YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC','WGQGTLVTVS']
pattern=re.compile('^'+'([ACDEFGHIKLMNPQRSTVWY]{5,30})'.join(map(re.escape,FRAMEWORK))+'$')
canonical=json.loads((ROOT/'P19438.json').read_text())['sequence']['value']
sources=[]
for name in ('bc2','bc2_v2','rfd4','rfantibody','rfantibody_crd3'):
    report=ROOT/'reports'/('canary_collection_'+name+'.json')
    if report.exists():
        sources.extend((name,Path(row['source_file'])) for row in json.loads(report.read_text()).get('candidates',[]))
sources.extend(('rfd4',p) for p in (ROOT/'generation/rfd4_canary/inference_samples').rglob('*_generated.cif'))
rows,failures,seen=[],[],set()
for method,path in sources:
    try:
        atoms=get_structure(CIFFile.read(path),model=1) if path.suffix=='.cif' else PDBFile.read(path).get_structure(model=1)
        atoms=atoms[np.isin(atoms.res_name,list(AA)) & (atoms.element!='H')]
        assert np.isfinite(atoms.coord).all()
        seqs={str(cid):''.join(AA[x] for x in atoms.res_name[(atoms.chain_id==cid)&(atoms.atom_name=='CA')]) for cid in np.unique(atoms.chain_id)}
        target_id=next(cid for cid,seq in seqs.items() if seq==canonical[29:211])
        binder_id=next(cid for cid in seqs if cid!=target_id)
        sequence=seqs[binder_id]
        key=(method,sequence)
        if key in seen: continue
        seen.add(key)
        match=pattern.fullmatch(sequence)
        assert match and 100<=len(sequence)<=150
        target=atoms[atoms.chain_id==target_id]
        binder=atoms[atoms.chain_id==binder_id]
        bca=binder[binder.atom_name=='CA']; tca=target[target.atom_name=='CA']
        d=np.linalg.norm(binder.coord[:,None,:]-target.coord[None,:,:],axis=-1)
        bmin=d.min(axis=1); tmin=d.min(axis=0)
        contact_b=np.isin(bca.res_id,np.unique(binder.res_id[bmin<=4.5]))
        cdrmask=np.zeros(len(sequence),dtype=bool)
        for group in (1,2,3): cdrmask[match.start(group):match.end(group)]=True
        cdr_share=float(cdrmask[contact_b].mean()) if contact_b.any() else None
        mutable=np.array([26<=i<=35 or 55<=i<=60 or 102<=i<=117 for i in range(1,len(sequence)+1)])
        mutable_share=float(mutable[contact_b].mean()) if method.startswith('bc2') and contact_b.any() else None
        anchor_clearance={}
        for pos in (54,145,151):
            rid=tca.res_id[pos-30]
            anchor=target.coord[(target.res_id==rid)&(target.atom_name=='ND2')]
            anchor_clearance[str(pos)]=float(np.linalg.norm(binder.coord-anchor[0],axis=1).min()) if len(anchor) else None
        backbone_only=bool(np.isin(binder.atom_name,['N','CA','C','O','OXT']).all())
        sasa=None if backbone_only else struc.sasa(binder,point_number=300)
        hydro=np.isin(binder.res_name,['ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP'])
        ca_dist=np.linalg.norm(np.diff(bca.coord,axis=0),axis=-1)
        target_positions=[30+i for i,rid in enumerate(tca.res_id) if rid in set(target.res_id[tmin<=4.5])]
        row={'method':method,'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),
             'source_file':str(path),'length':len(sequence),'CDR_sequences':list(match.groups()),
             'interface_binder_residue_count_4p5A':int(contact_b.sum()),
             'interface_target_uniprot_positions_4p5A':target_positions,
             'interface_full_CDR_residue_fraction_4p5A':cdr_share,
             'interface_BC2_editable_paratope_fraction_4p5A':mutable_share,
             'glycan_ND2_minimum_binder_distance_A':anchor_clearance,
             'source_atom_scope':'backbone_only' if backbone_only else 'all_or_partial_sidechains',
             'free_binder_hydrophobic_SASA_fraction':None if backbone_only else float(np.nansum(sasa[hydro])/np.nansum(sasa)),
             'consecutive_CA_distance_min_A':float(ca_dist.min()),'consecutive_CA_distance_max_A':float(ca_dist.max()),
             'cysteine_count':sequence.count('C'),
             'longest_homopolymer':max(len(list(g)) for _,g in itertools.groupby(sequence)),
             'potential_N_glycosylation_sequons':[{'position':m.start()+1,'motif':m.group()} for m in re.finditer('N[^P][ST]',sequence)],
             'root_backbone_name':path.stem.split('_candidate')[0] if method.startswith('bc2') else path.stem}
        rows.append(row)
    except Exception as exc: failures.append({'source_file':str(path),'error':repr(exc)})
summary={'job_id':os.environ['SLURM_JOB_ID'],'rows':rows,'failures':failures,
         'scope':'CPU structural descriptors; full CDR segments use the 83-residue source-core definition. Glycan-anchor clearance does not model complete glycans. SASA is not a solubility measurement.'}
(ROOT/'reports/parallel_quality_audit_v2.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({'audited_sequences':len(rows),'failures':len(failures),
                  'methods':sorted({r['method'] for r in rows})}),flush=True)
if not rows: raise SystemExit('No cases audited')
