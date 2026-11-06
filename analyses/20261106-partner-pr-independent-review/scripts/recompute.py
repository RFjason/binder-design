"""Independent PR checks against frozen results and public per-seed data."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import statistics as st
import subprocess
import sys
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--local-only',action='store_true');parser.add_argument('--data-dir');parser.add_argument('--dependency-dir');args=parser.parse_args()
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parents[1];SOURCES=PACKAGE/'sources'
ROOT=PACKAGE/'artifacts/recomputed';ROOT.mkdir(parents=True,exist_ok=True)
DATA=Path(args.data_dir).resolve() if args.data_dir else PACKAGE/'external_data'
if args.dependency_dir:sys.path.insert(0,str(Path(args.dependency_dir).resolve()))
rows=[json.loads(s) for s in (REPO/'campaigns/tnfr1-vhh-20261103/results/all_candidates.jsonl').read_text().splitlines()]
spec=importlib.util.spec_from_file_location('peer_check',SOURCES/'pr1/tools/check_submission.py');tool=importlib.util.module_from_spec(spec);spec.loader.exec_module(tool)
hash_checks=[]
for line in (SOURCES/'pr1/analyses/20261105-compliance-novelty-audit/input_hashes.txt').read_text().splitlines():
 digest,name=line.split(None,1);p=(SOURCES/'pr1'/name) if name.startswith('analyses/') else (REPO/name)
 hash_checks.append({'path':name,'matches':hashlib.sha256(p.read_bytes()).hexdigest()==digest})
assert all(v['matches'] for v in hash_checks)
route=[];patch=[]
for method in ('rfd4','rfantibody','bindcraft2'):
 group=[r for r in rows if r['structure_method']==method];hits=[r for r in group if r['computational_quality_pass']]
 route.append({'method':method,'sequences':len(group),'source_parents_with_output':len({r['parent_id'] for r in group}),'passing_sequences':len(hits),'passing_parents':len({r['parent_id'] for r in hits}),
               'ipsae_ge_0p6':sum(r['ipsae_min']>=.6 for r in group),'ipsae_pass_but_pose_fail':sum(r['ipsae_min']>=.6 and (r.get('pose_dockq') is None or r['pose_dockq']<.23) for r in group),'pose_dockq_median':st.median([r['pose_dockq'] for r in group if r.get('pose_dockq') is not None])})
 for p in ('crd2','crd3'):
  g=[r for r in group if r['stratum'].startswith(p)];h=[r for r in g if r['computational_quality_pass']]
  patch.append({'method':method,'patch':p,'sequences':len(g),'parents_with_output':len({r['parent_id'] for r in g}),'passing_sequences':len(h),'passing_parents':len({r['parent_id'] for r in h})})
# This regression case isolates the omission of CDR2 from the key; numbering
# is supplied explicitly so dependency availability cannot alter the result.
seq1='QVQLVESGGGLVQPGGSLRLSCAASGDISNITYWGWFRQAPGQGLEAVAYISPGGDINYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAAGSSGAFDYWGQGTLVTVSS'
seq2=seq1.replace('YISPGGDIN','YISPGADIN')
tool._CDR_CACHE[seq1]={'CDR1':'GDISNITY','CDR2':'YISPGGDIN','CDR3':'AAGSSGAFDY'}
tool._CDR_CACHE[seq2]={'CDR1':'GDISNITY','CDR2':'YISPGADIN','CDR3':'AAGSSGAFDY'}
false_pairs=tool.dedup_within([('q1',seq1),('q2',seq2)],'vhh');assert false_pairs
# An unavailable requested novelty backend is currently silently treated as
# zero hits, leaving the row PASS and process exit0.
fasta=ROOT/'novelty_backend_repro.fasta';fasta.write_text('>q1\n'+seq1+'\n')
ref=SOURCES/'pr1/analyses/20261105-compliance-novelty-audit/reference_library/TNFR1_known_binders.fasta'
call=subprocess.run([sys.executable,str(SOURCES/'pr1/tools/check_submission.py'),str(fasta),'vhh','--no-anarcii','--refs',str(ref),'--mmseqs',str(ROOT/'definitely_absent_mmseqs'),'--out',str(ROOT/'novelty_backend_repro.csv')],capture_output=True,text=True,encoding='utf8',env={**__import__('os').environ,'PYTHONIOENCODING':'utf-8'})
assert call.returncode==0
(ROOT/'novelty_backend_repro_stdout.txt').write_text(call.stdout.replace(str(ROOT),'[OUTPUT_DIR]'),encoding='utf8')
report={'input_hash_checks':hash_checks,'route':route,'patch':patch,'cdr2_omission_repro':false_pairs,'missing_novelty_backend_exit_code':call.returncode}
(ROOT/'local_checks.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:v for k,v in report.items() if k!='input_hash_checks'},indent=2))

# Reproduce the external table using a declared seed0 and design stoichiometry.
if args.local_only:raise SystemExit(0)
import pyarrow.parquet as pq
d=pq.read_table(DATA/'cofold_predictions.parquet').to_pandas();s=pq.read_table(DATA/'design_summary.parquet').to_pandas()
d=d[(d['cofolding_model']=='ef2full') & d['target_form'].fillna('').eq('')].copy()
desired={name:'1to2' for name in ('15-PGDH','VEGF-A','Mature GDF-8','Latent GDF-8')};desired.update({'TNFa':'1to3','Cas9':'rnp'})
d=d[d.apply(lambda r:r.stoichiometry==desired.get(r.target,'1to1'),axis=1)]
j=d.merge(s[['uuid','binder_final']],on='uuid',validate='many_to_one');j=j[j.binder_final.notna()]
one=j[j.seed.astype(str)=='0'].copy();assert len(one)==1320 and one.uuid.nunique()==1320
max_scores=j.groupby('uuid').ipsae_min.max();delta=one.uuid.map(max_scores)-one.ipsae_min
calibration=[]
for threshold in (.5,.55,.6,.65,.7,.8):
 passed=one.ipsae_min>=threshold;pos=one.binder_final.astype(bool)
 calibration.append({'threshold':threshold,'TP':int((passed&pos).sum()),'FP':int((passed&~pos).sum()),'FN':int((~passed&pos).sum()),'TN':int((~passed&~pos).sum())})
single=one.set_index('uuid').ipsae_min
recheck={'seed':0,'n':len(one),'positive':int(one.binder_final.sum()),'negative':int((~one.binder_final.astype(bool)).sum()),'positive_median':float(one[one.binder_final.astype(bool)].ipsae_min.median()),'negative_median':float(one[~one.binder_final.astype(bool)].ipsae_min.median()),'five_seed_max_minus_seed0_median':float(delta.median()),'five_seed_max_minus_seed0_mean':float(delta.mean()),'five_seed_gap_gt_0p1':int((delta>.1).sum()),'five_seed_changes_ge_0p6':int(((single>=.6)!=(max_scores.reindex(single.index)>=.6)).sum()),'calibration':calibration,'source_sha256':{n:hashlib.sha256((DATA/n).read_bytes()).hexdigest() for n in ('cofold_predictions.parquet','design_summary.parquet')}}
(ROOT/'calibration_recomputed.json').write_text(json.dumps(recheck,indent=2));print('CALIBRATION',json.dumps(recheck))
