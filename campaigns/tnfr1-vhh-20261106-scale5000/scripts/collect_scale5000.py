"""One serialized collector owns quota registration and exact-sequence de-duplication."""
import hashlib
import json
import os
from pathlib import Path
import re
import traceback

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile,get_structure
from s5k_common import atomic_json,read_json,registry_lock,select_row,progress

assert os.environ.get('SLURM_JOB_ID'),'Slurm only'
CONFIG=Path(__file__).resolve().parents[1]/'config/campaign.json'
campaign=read_json(CONFIG);ROOT=Path(campaign['remote_root'])
b=campaign['batches'][int(os.environ['S5K_TASK'])];label='scale5000_'+b['batch_id'];source=Path(b['output_dir'])
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
core=['QVQLVESGGGLVQPGGSLRLSC','WFRQAPGQGLEAVA','YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC','WGQGTLVTVS']
contract=read_json(b['framework_contract']) if b.get('framework_contract') else None
if contract:core=contract['core']
pattern=re.compile('([ACDEFGHIKLMNPQRSTVWY]{5,30})'.join(map(re.escape,core)))
canonical=read_json(ROOT/'P19438.json')['sequence']['value']
if b['method']=='rfd4':files=sorted((source/'inference_samples').rglob('*_generated.cif'))
elif b['method']=='bc2_v2':
    files=sorted(p for p in source.rglob('*.cif') if 'Complexes' in p.parts)
    files+=sorted((source/'1_Trajectories').rglob('trajectory.cif'))
else:files=sorted((source/'design_pdbs').glob('*.pdb'))
rows=[];rejects=[];duplicates=[];quota_skipped=[]
with registry_lock(ROOT):
    registry=read_json(ROOT/'state/registry.json',{})
    for path in files:
        owned_key=None
        try:
            atoms=get_structure(CIFFile.read(path),model=1) if path.suffix=='.cif' else PDBFile.read(path).get_structure(model=1)
            atoms=atoms[np.isin(atoms.res_name,list(AA)) & (atoms.element!='H')]
            assert np.isfinite(atoms.coord).all(),'nonfinite coordinates'
            seqs={str(cid):''.join(AA[x] for x in atoms.res_name[(atoms.chain_id==cid)&(atoms.atom_name=='CA')]) for cid in np.unique(atoms.chain_id)}
            target=[cid for cid,seq in seqs.items() if seq==canonical[29:211]];assert len(target)==1
            target=target[0];binder_ids=[cid for cid in seqs if cid!=target];assert len(binder_ids)==1
            binder_id=binder_ids[0];seq=seqs[binder_id];match=pattern.fullmatch(seq)
            assert match and seq.count('C')%2==0
            source_sequence=seq
            # The coordinate scaffold ends at IMGT127. Score the complete mature
            # sequence including S128; do not invent coordinates in the reference.
            seq=seq if contract else seq+'S'
            assert 100<=len(seq)<=150
            digest=hashlib.sha256(seq.encode()).hexdigest();cid='TNFR1_s5k_'+b['batch_id']+'_'+digest[:12]
            binder=atoms[atoms.chain_id==binder_id];ca=binder[binder.atom_name=='CA']
            sg=binder[(binder.res_name=='CYS')&(binder.atom_name=='SG')]
            geometry=None;core_distance=None
            if len(sg)==seq.count('C'):
                dist=np.linalg.norm(sg.coord[:,None,:]-sg.coord[None,:,:],axis=-1);np.fill_diagonal(dist,np.inf)
                res=[ca.res_id[contract['cys1_index0'] if contract else len(core[0])-1],ca.res_id[match.start(3)-1]]
                indices=[list(sg.res_id).index(x) for x in res];core_distance=float(dist[indices[0],indices[1]])
                def pairable(rest):
                    if not rest:return True
                    return any(dist[rest[0],p]<=3 and pairable([v for v in rest[1:] if v!=p]) for p in rest[1:])
                geometry=bool(core_distance<=3 and pairable([i for i in range(len(sg)) if i not in indices]))
                assert geometry,'source disulfide geometry failed'
            parent=(b['parent_id'] if b['method']=='mpnn_refine' else b['batch_id']+':trajectory_1' if b['method']=='bc2_v2' else b['batch_id']+':'+path.stem.split('_dldesign_')[0])
            row={'design_id':cid,'sequence':seq,'length':len(seq),'method':b['method'],
                 'structure_method':b.get('parent_structure_method',{'rfd4':'rfd4','bc2_v2':'bindcraft2','rfantibody':'rfantibody'}.get(b['method'])),
                 'seq_method':'proteinmpnn' if b['method'] in ('rfantibody','mpnn_refine') else b['method'],
                 'campaign':'scale5000','batch_id':b['batch_id'],'stratum':b['stratum'],'parent_id':parent,
                 'optimization_parent':b.get('parent_design_id'),'generation_seed':b['seed'],
                 'source_file':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                 'canonical_target_verified':True,'framework_identity_verified':True,'frame_endpoints_checked':True,
                 'fixed_framework_residues':sum(map(len,core)) if contract else 84,'cdr_sequences':list(match.groups()),'cdr_lengths':list(map(len,match.groups())),
                 'framework_id':b.get('framework_id','h-NbBCII10'),
                 'framework_core':core if contract else core[:-1]+[core[-1]+'S'],
                 'framework_cys1_index0':contract['cys1_index0'] if contract else len(core[0])-1,
                 'cdr_definition':'IMGT' if contract else 'legacy source design segments',
                 'source_sequence':source_sequence,'maturation':'complete_IMGT_domain_before_design' if contract else 'append_IMGT_S128_before_folding',
                 'reference_terminal_S128_coordinates':'present' if contract else 'absent; DockQ uses aligned resolved residues',
                 'cysteine_count':seq.count('C'),'cysteine_geometry_ok':geometry,'core_cysteine_SG_distance_A':core_distance,
                 'source_atom_scope':'backbone_only' if set(binder.atom_name)<=set(['N','CA','C','O','CB']) else 'all_or_partial_sidechains',
                 'source_SG_atom_count':len(sg),'novelty_status':'NOT_RUN','status':'PENDING_FROZEN_SCORING',
                 'final_rank_eligible':False}
            decision=select_row(registry,row,b['batch_id'],campaign['new_quota'])
            if decision=='duplicate':duplicates.append({'source':str(path),'sequence_sha256':digest});continue
            if decision=='quota_full':quota_skipped.append(str(path));continue
            if decision=='new':owned_key=digest
            folder=ROOT/'candidates'/cid;folder.mkdir(parents=True,exist_ok=True)
            reference=atoms.copy();original=reference.chain_id.copy();reference.chain_id[original==target]='T';reference.chain_id[original==binder_id]='H'
            reference.set_annotation('b_factor',np.zeros(len(reference),dtype=float))
            pdb=PDBFile();pdb.set_structure(reference);pdb.write(folder/'designed_complex.pdb')
            row['reference_pdb']=str(folder/'designed_complex.pdb')
            (folder/'sequence.fasta').write_text('>'+cid+'\n'+seq+'\n');atomic_json(folder/'provenance.json',row)
            if not any(r['design_id']==cid for r in rows):rows.append(row)
        except Exception as exc:
            if owned_key is not None:registry.pop(owned_key,None)
            rejects.append({'source':str(path),'reason':repr(exc),'traceback':traceback.format_exc()})
    cases=[]
    for row in rows:
        cases.extend([{'id':row['design_id']+'_complex','candidate_id':row['design_id'],'chains':[{'id':'A','sequence':canonical[:211],'msa':'inputs/TNFR1_competition_M1_T211.a3m'},{'id':'B','sequence':row['sequence']}]},
                      {'id':row['design_id']+'_monomer','candidate_id':row['design_id'],'chains':[{'id':'B','sequence':row['sequence']}]}])
    manifest=ROOT/'inputs'/('canary_scoring_'+label+'.json');atomic_json(manifest,cases)
    atomic_json(ROOT/'reports'/('canary_collection_'+label+'.json'),{'method':b['method'],'batch':b,'source_files':len(files),'candidates':rows,'rejected_files':rejects,'duplicates_skipped':duplicates,'quota_skipped':quota_skipped,'scoring_manifest':str(manifest)})
    atomic_json(ROOT/'state/registry.json',registry)
progress(ROOT,campaign)
print(json.dumps({'batch_id':b['batch_id'],'collected':len(rows),'rejected':len(rejects),'duplicates':len(duplicates),'quota_skipped':len(quota_skipped)}),flush=True)
