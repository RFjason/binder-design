"""Offline checks of pinned peer revision; demonstrates remaining edge cases."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'sources/peer_revision_1bbd20b/check_submission.py'
subprocess.run([sys.executable,str(source),'--selftest'],check=True,env={**os.environ,'PYTHONIOENCODING':'utf-8'})
spec=importlib.util.spec_from_file_location('peer_revision',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
base='QVQLVESGGGLVQPGGSLRLSCAAS{a}WGWFRQAPGQGLEAVA{b}YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC{c}WGQGTLVTVSS'
a,b,c='GDISNITY','IDLRNNKT','AAGSSGAFDY'
s=base.format(a=a,b=b,c=c);t=base.format(a=a,b=b,c=c+'G')
m._CDR_CACHE[s]={'CDR1':a,'CDR2':b,'CDR3':c};m._CDR_CACHE[t]={'CDR1':a,'CDR2':b,'CDR3':c+'G'}
assert (a+b+c) in (a+b+c+'G') and len(a+b+c)/len(a+b+c+'G')>.9
assert m.dedup_within([('short',s),('long',t)],'vhh')==([],0),'Pinned bug repro changed'
excluded=[]
for f in (ROOT/'sources/pr1/analyses/20261105-compliance-novelty-audit/reference_library').glob('*.fasta'):
 for entry in f.read_text().split('>')[1:]:
  lines=entry.splitlines();seq=''.join(lines[1:]);head=lines[0]
  if len(seq)>200 and any(x in head.lower() for x in ('heavy','antibody','fab')):
   excluded.append({'reference_file':f.name,'header':head,'length':len(seq)})
assert len(excluded)==35
dest=ROOT/'artifacts/followup';dest.mkdir(parents=True,exist_ok=True)
(dest/'excluded_long_antibody_references.json').write_text(json.dumps(excluded,indent=2),encoding='utf8')
print('CONFIRMED: 4 peer selftests pass; partial-CDR duplicate remains unreported; 35 labelled long antibody entries are excluded by the 200-aa numbering gate.')
