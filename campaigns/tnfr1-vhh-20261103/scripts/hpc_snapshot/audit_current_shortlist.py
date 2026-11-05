"""Quantitative monomer/CDR descriptors for all current fixed-protocol screening hits."""
import json
import os
from pathlib import Path
import re
import time

import numpy as np
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile

assert os.environ.get('SLURM_JOB_ID'), 'Compute node only'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
core=['QVQLVESGGGLVQPGGSLRLSC','WFRQAPGQGLEAVA','YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC','WGQGTLVTVS']
pattern=re.compile('([ACDEFGHIKLMNPQRSTVWY]{5,30})'.join(map(re.escape,core)))
candidates={}
for glob in ('canary_result_bulk1_*.json','canary_result_bc2_recovery_*.json'):
    for p in (ROOT/'reports').glob(glob):
        try:report=json.loads(p.read_text())
        except ValueError:continue
        for row in report['scored_candidates']:
            if row['computational_quality_pass']:candidates[row['design_id']]=row
rows,failures=[],[]
for row in candidates.values():
    try:
        seq=row['sequence'];match=pattern.fullmatch(seq);assert match
        file=ROOT/'predictions'/(row['design_id']+'_monomer')/'model.pdb'
        atoms=PDBFile.read(file).get_structure(model=1,extra_fields=['b_factor'])
        atoms=atoms[(atoms.chain_id=='B') & (atoms.element!='H')]
        ca=atoms[atoms.atom_name=='CA'];assert len(ca)==len(seq)
        confidence=ca.b_factor.astype(float)
        if confidence.max()<=1.1:confidence*=100
        sasa=struc.sasa(atoms,point_number=600)
        assert np.isfinite(sasa).any()
        total_sasa=float(np.nansum(sasa));assert total_sasa>0
        hydrophobic=np.isin(atoms.res_name,['ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP'])
        cdrs=[]
        for i in (1,2,3):
            start,end=match.span(i);resids=ca.res_id[start:end]
            part=np.isin(atoms.res_id,resids);part_sasa=float(np.nansum(sasa[part]))
            cdrs.append({'CDR':i,'sequence':seq[start:end],'start_1based':start+1,'end_1based':end,
                         'mean_plddt':float(confidence[start:end].mean()),
                         'minimum_plddt':float(confidence[start:end].min()),
                         'fraction_plddt_below_70':float((confidence[start:end]<70).mean()),
                         'SASA_A2':part_sasa,'hydrophobic_SASA_A2':float(np.nansum(sasa[part & hydrophobic]))})
        sequons=[]
        for m in re.finditer('(?=(N[^P][ST]))',seq):
            pos=m.start();sequons.append({'position_1based':pos+1,'motif':m.group(1),
                                         'CDR':next((i for i in (1,2,3) if match.start(i)<=pos<match.end(i)),None)})
        rows.append({**row,'monomer_total_SASA_A2':total_sasa,
                     'monomer_hydrophobic_SASA_fraction':float(np.nansum(sasa[hydrophobic])/total_sasa),
                     'CDR_local_descriptors':cdrs,'potential_N_glycosylation_sequons':sequons,
                     'longest_homopolymer':max(len(m.group()) for m in re.finditer(r'(.)\1*',seq)),
                     'novelty_audit':'NOT_RUN','full_glycan_accessibility_audit':'NOT_RUN',
                     'final_rank_eligible':False})
    except Exception as exc:failures.append({'design_id':row['design_id'],'error':repr(exc)})
rows.sort(key=lambda r:r['ipsae_min'],reverse=True)
out={'job_id':os.environ['SLURM_JOB_ID'],'updated_epoch':time.time(),'candidate_count':len(rows),
     'distinct_parent_ids':len({r['parent_id'] for r in rows}),'candidates':rows,'failures':failures,
     'scope':'Sequence liabilities and predicted monomer surface/CDR confidence descriptors. Final selection still requires novelty and full glycan accessibility audits.'}
(ROOT/'reports/current_shortlist_descriptors.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'candidate_count':len(rows),'distinct_parent_ids':out['distinct_parent_ids'],'failures':failures}),flush=True)
if failures:raise SystemExit(1)
