"""Package immutable campaign results for the owner's private collaboration repository."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
import zipfile

assert os.environ.get('SLURM_JOB_ID'), 'Export I/O and compression must use a compute node'
ROOT=Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
EXPORT=ROOT/'exports/github_tnfr1_v1'
EXPORT.mkdir(parents=True,exist_ok=True)
CURATED=EXPORT/'curated';CURATED.mkdir(exist_ok=True)
def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def copy(source,dest):
    source=Path(source);dest=CURATED/dest
    assert source.is_file(),str(source)
    dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
    return dest
def write_json(name,value):
    p=CURATED/name;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')

rows={};source_reports=[]
for pattern in ('canary_result_bulk1_*.json','canary_result_bc2_recovery_*.json'):
    for p in sorted((ROOT/'reports').glob(pattern)):
        report=json.loads(p.read_text());source_reports.append({'name':p.name,'sha256':digest(p)})
        for row in report['scored_candidates']:
            cid=row['design_id'];assert cid not in rows,cid
            rows[cid]={**row,'source_report':p.name}
assert len(rows)==1127,len(rows)
ordered=sorted(rows.values(),key=lambda r:(not r['computational_quality_pass'],-r['ipsae_min'],r['design_id']))
hits=[r for r in ordered if r['computational_quality_pass']]
assert len(hits)==19 and len({r['parent_id'] for r in hits})==8
protocol=json.loads((ROOT/'scoring_protocol_v1.json').read_text())
all_artifacts=[]
for i,row in enumerate(ordered):
    cid=row['design_id'];seq=row['sequence']
    assert 100<=len(seq)<=150 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY')
    assert hashlib.sha256(seq.encode()).hexdigest()[:12]==cid.rsplit('_',1)[-1]
    assert not row['final_rank_eligible']
    for suffix in ('complex','monomer'):
        folder=ROOT/'predictions'/(cid+'_'+suffix)
        s=json.loads((folder/'summary.json').read_text())
        assert s['ok'] and s['kwargs']['seed']==0 and s['kwargs']['num_diffusion_samples']==1
        assert s['templates'] is False and s['initial_guess'] is False
        if suffix=='complex':
            assert s['protocol_id']==protocol['protocol_id']
            assert s['msa_records'][0]['sha256']==protocol['target']['msa_sha256']
            assert s['msa_records'][1]['rows']==0
            assert s['ipsae_min']==row['ipsae_min']
        for filename in ('model.cif','pae.npz','summary.json','input_and_token_mapping.json'):
            path=folder/filename;assert path.is_file(),str(path)
            all_artifacts.append((path,f'predictions/{cid}/{suffix}/{filename}'))
    all_artifacts.append((Path(row['reference_pdb']),f'designs/{cid}/designed_complex.pdb'))
    original=Path(row['source_file'])
    assert digest(original)==row['source_sha256']
    all_artifacts.append((original,f'designs/{cid}/original{original.suffix}'))
    if i%100==0:print('EXPORT_VALIDATED',i,flush=True)

(CURATED/'results').mkdir(exist_ok=True)
with (CURATED/'results/all_candidates.jsonl').open('w') as f:
    for row in ordered:f.write(json.dumps(row,ensure_ascii=False)+'\n')
columns=['design_id','structure_method','campaign','stratum','parent_id','length','sequence',
         'cdr1','cdr2','cdr3','ipsae_min','ipTM','monomer_plddt','pose_dockq','pose_dockq_status',
         'monomer_disulfide_geometry_ok','predicted_full_CDR_interface_fraction_4p5A',
         'predicted_signal_peptide_contact_binder_residues_4p5A','computational_quality_pass',
         'final_rank_eligible','novelty_status','source_sha256','source_report']
def table(name,items):
    path=CURATED/name;path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');writer.writeheader()
        for row in items:
            writer.writerow({**row,**dict(zip(('cdr1','cdr2','cdr3'),row['cdr_sequences']))})
table('results/all_candidates.csv',ordered);table('shortlist/candidates.csv',hits)
(CURATED/'shortlist/candidates.fasta').write_text(''.join('>'+r['design_id']+' COMPUTATIONAL_SCREEN_ONLY\n'+r['sequence']+'\n' for r in hits))
for row in hits:
    cid=row['design_id']
    copy(row['reference_pdb'],f'shortlist/structures/{cid}/designed_complex.pdb')
    for suffix in ('complex','monomer'):
        folder=ROOT/'predictions'/(cid+'_'+suffix)
        copy(folder/'model.cif',f'shortlist/structures/{cid}/{suffix}.cif')
        copy(folder/'summary.json',f'shortlist/structures/{cid}/{suffix}_summary.json')

for name in ('TNFR1_competition_M1_T211.fasta','TNFR1_mature_L30_T211.fasta',
             'TNFR1_mature_canonical.pdb','TNFR1_competition_M1_T211.a3m',
             'h-NbBCII10_complete.cif','hotspots_crd2.txt','hotspots_crd3.txt'):
    copy(ROOT/'inputs'/name,'inputs/'+name)
for p in (ROOT/'inputs/rfd4').glob('*.cif'):copy(p,'inputs/rfd4/'+p.name)
for name in ('bulk1','bc2_recovery'):
    for p in (ROOT/'inputs'/name).glob('*.json'):copy(p,f'inputs/{name}/{p.name}')
copy(ROOT/'scoring_protocol_v1.json','config/scoring_protocol_v1.json')
copy(ROOT/'P19438.json','inputs/uniprot_P19438.json')
for p in (ROOT/'scripts').iterdir():
    if p.suffix in ('.py','.sh','.slurm'):copy(p,'scripts/hpc_snapshot/'+p.name)
missing_optional=[]
for name in ('job_ledger.json','target_planning_audit.json','local_structure_audit.json',
             'design_input_audit.json','bulk1_preflight.json','bulk1_progress.json',
             'bc2_recovery_progress.json','bc2_gpu_backend_probe.json','current_shortlist_descriptors.json'):
    if (ROOT/'reports'/name).is_file():copy(ROOT/'reports'/name,'provenance/'+name)
    else:missing_optional.append(name)
write_json('provenance/hpc_export_notes.json',{'optional_reports_not_present_on_hpc':missing_optional,
    'note':'Planning reports made on the Windows workspace may be added separately with their own hashes.'})
write_json('provenance/source_reports.json',source_reports)
methods={name:{'scored':sum(r['structure_method']==name for r in ordered),
               'screen_pass':sum(r['structure_method']==name for r in hits)}
         for name in ('rfd4','rfantibody','bindcraft2')}
write_json('results/summary.json',{'campaign':'tnfr1-vhh-20261103','export_job':os.environ['SLURM_JOB_ID'],
    'export_epoch':time.time(),'scored_candidates':len(rows),'screen_pass_candidates':len(hits),
    'passing_parent_trajectories':8,'methods':methods,'final_submission_eligible':0,
    'status':'COMPUTATIONAL_SCREENING_ONLY','allocated_gpu_hours':74.19666666666662,
    'pending_final_audits':['public_sequence_novelty','full_target_glycan_accessibility','diversity_and_final_selection']})

# Split raw prediction outputs by uncompressed bytes; no model weights or credentials.
max_part_bytes=600*1024*1024
parts=[];archive=None;part_bytes=0
def close_part():
    if archive is not None:
        archive.close();parts[-1].update(bytes=part_path.stat().st_size,sha256=digest(part_path))
try:
    for source,relative in all_artifacts:
        size=source.stat().st_size
        if archive is None or part_bytes+size>max_part_bytes:
            close_part();part_path=EXPORT/f'tnfr1-full-predictions-part{len(parts)+1:02d}.zip'
            archive=zipfile.ZipFile(part_path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3)
            parts.append({'filename':part_path.name,'files':0,'uncompressed_bytes':0});part_bytes=0
            print('EXPORT_ARCHIVE',part_path.name,flush=True)
        archive.write(source,relative);part_bytes+=size
        parts[-1]['files']+=1;parts[-1]['uncompressed_bytes']+=size
finally:close_part()
write_json('provenance/archive_manifest.json',{'status':'FULL_DESIGN_AND_PREDICTION_OUTPUTS',
    'candidate_count':1127,'files':len(all_artifacts),'parts':parts,
    'included_per_candidate':['original_generated_structure','designed_complex.pdb',
       'complex/model.cif','complex/pae.npz','complex/summary.json','complex/input_and_token_mapping.json',
       'monomer/model.cif','monomer/pae.npz','monomer/summary.json','monomer/input_and_token_mapping.json']})
files={str(p.relative_to(CURATED)):{'bytes':p.stat().st_size,'sha256':digest(p)}
       for p in sorted(CURATED.rglob('*')) if p.is_file() and p.name!='data_manifest.json'}
write_json('provenance/data_manifest.json',{'files':files,'file_count':len(files)})
bundle=EXPORT/'tnfr1-collaboration-bundle.zip'
with zipfile.ZipFile(bundle,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3) as z:
    for p in sorted(CURATED.rglob('*')):
        if p.is_file():z.write(p,str(p.relative_to(CURATED)))
status={'ok':True,'job_id':os.environ['SLURM_JOB_ID'],'export_root':str(EXPORT),
        'curated_bundle':{'filename':bundle.name,'bytes':bundle.stat().st_size,'sha256':digest(bundle)},
        'raw_parts':parts,'candidate_count':1127,'shortlist_count':19}
(EXPORT/'export_status.json').write_text(json.dumps(status,indent=2)+'\n')
print('EXPORT_COMPLETE',json.dumps(status),flush=True)
