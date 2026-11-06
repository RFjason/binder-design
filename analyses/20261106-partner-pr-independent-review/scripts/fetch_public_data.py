"""Fetch the exact public data revision used here and require its SHA-256."""
import hashlib
import json
from pathlib import Path
import urllib.request
ROOT=Path(__file__).resolve().parents[1]
meta=json.loads((ROOT/'manifest.json').read_text(encoding='utf8'))
dest=ROOT/'external_data';dest.mkdir(exist_ok=True)
for name,expected in meta['public_data_sha256'].items():
 path=dest/name
 if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==expected:
  print('ALREADY_VERIFIED',name);continue
 folder='insilico/' if name=='cofold_predictions.parquet' else ''
 url='https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/resolve/'+meta['public_dataset_revision']+'/data/tables/'+folder+name
 temp=path.with_suffix('.download')
 urllib.request.urlretrieve(url,temp)
 actual=hashlib.sha256(temp.read_bytes()).hexdigest()
 if actual!=expected:raise RuntimeError(f'Public input hash mismatch for {name}: {actual} != {expected}')
 temp.replace(path);print('VERIFIED',name,expected)
