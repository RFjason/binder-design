"""Small shared state primitives; no model imports."""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import time
import uuid

def read_json(path,default=None):
    try:return json.loads(Path(path).read_text())
    except (FileNotFoundError,json.JSONDecodeError):return default

def atomic_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    temp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n');temp.replace(path)

@contextmanager
def registry_lock(root):
    path=Path(root)/'state/registry.lock';path.parent.mkdir(parents=True,exist_ok=True)
    deadline=time.monotonic()+60
    while True:
        try:path.mkdir();break
        except FileExistsError:
            if time.monotonic()>deadline:raise RuntimeError('Registry lock remains held; inspect its owner before recovery')
            time.sleep(0.2)
    try:
        atomic_json(path/'owner.json',{'job_id':os.environ.get('SLURM_JOB_ID'),'pid':os.getpid(),'time':time.time()})
        yield
    finally:
        (path/'owner.json').unlink(missing_ok=True);path.rmdir()

def new_count(registry):return sum(v['origin']=='new' for v in registry.values())

def select_row(registry,row,batch_id,quota):
    key=hashlib.sha256(row['sequence'].encode()).hexdigest()
    existing=registry.get(key)
    if existing:
        return 'existing_own' if existing.get('design_id')==row['design_id'] else 'duplicate'
    if new_count(registry)>=quota:return 'quota_full'
    registry[key]={'origin':'new','design_id':row['design_id'],'batch_id':batch_id,'parent_id':row['parent_id']}
    return 'new'

def progress(root,campaign):
    root=Path(root);registry=read_json(root/'state/registry.json',{})
    scored={};failures=[];collected_batches=0
    for batch in campaign['batches']:
        label='scale5000_'+batch['batch_id']
        c=read_json(root/'reports'/('canary_collection_'+label+'.json'))
        if c is not None:collected_batches+=1
        report=read_json(root/'reports'/('canary_result_'+label+'.json'))
        if report:
            scored.update({r['design_id']:r for r in report['scored_candidates']})
            failures.extend(report['failures'])
    value={'campaign':'scale5000','updated_epoch':time.time(),'baseline_scored':1127,
           'new_reserved':new_count(registry),'new_scored_and_audited':len(scored),
           'cumulative_scored':1127+len(scored),'target':5000,'collection_batches':collected_batches,
           'new_quality_pass':sum(r['computational_quality_pass'] for r in scored.values()),
           'baseline_quality_pass':19,'audit_failures':failures,
           'final_submission_eligible':0}
    atomic_json(root/'reports/progress.json',value);return value
