"""Standard-library integrity checks; no protein inference or third-party dependencies."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

REPO=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--campaign',type=Path,default=REPO/'campaigns/tnfr1-vhh-20261103')
parser.add_argument('--archives',type=Path)
args=parser.parse_args();root=args.campaign.resolve()
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
manifest_paths=[root/'provenance/data_manifest.json']
if (root/'provenance/workspace_manifest.json').exists():manifest_paths.append(root/'provenance/workspace_manifest.json')
verified_files=0
for manifest_path in manifest_paths:
    manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    for name,expected in manifest['files'].items():
        path=(root/name).resolve();assert path.is_relative_to(root),name
        assert path.is_file(),f'Missing: {name}'
        assert path.stat().st_size==expected['bytes'],f'Byte count: {name}'
        assert sha(path)==expected['sha256'],f'SHA256: {name}'
        verified_files+=1
rows=[json.loads(line) for line in (root/'results/all_candidates.jsonl').read_text(encoding='utf-8').splitlines() if line]
assert len(rows)==len({r['design_id'] for r in rows})==1127
framework=['QVQLVESGGGLVQPGGSLRLSC','WFRQAPGQGLEAVA','YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC','WGQGTLVTVS']
pattern=re.compile('([ACDEFGHIKLMNPQRSTVWY]{5,30})'.join(map(re.escape,framework)))
for row in rows:
    seq=row['sequence'];assert 100<=len(seq)<=150 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY')
    assert row['length']==len(seq) and pattern.fullmatch(seq),row['design_id']
    assert hashlib.sha256(seq.encode()).hexdigest()[:12]==row['design_id'].rsplit('_',1)[-1]
    assert not row['final_rank_eligible']
    assert abs(row['ipsae_min']-min(row['ipsae_directions'].values()))<1e-6
    passing=bool(row['monomer_plddt']>=70 and row['pose_dockq'] is not None and row['pose_dockq']>=0.23
                 and row['monomer_disulfide_geometry_ok'] and row['ipsae_min']>=0.6)
    assert row['computational_quality_pass']==passing,row['design_id']
hits=[r for r in rows if r['computational_quality_pass']]
assert len(hits)==19 and len({r['parent_id'] for r in hits})==8
assert Counter(r['structure_method'] for r in rows)=={'rfd4':520,'rfantibody':357,'bindcraft2':250}
assert Counter(r['structure_method'] for r in hits)=={'rfd4':4,'bindcraft2':15}
fasta={};key=None
for line in (root/'shortlist/candidates.fasta').read_text(encoding='utf-8').splitlines():
    if line.startswith('>'):key=line[1:].split()[0];assert key not in fasta;fasta[key]=''
    elif line.strip():assert key is not None;fasta[key]+=line.strip()
assert fasta=={r['design_id']:r['sequence'] for r in hits}
protocol=json.loads((root/'config/scoring_protocol_v1.json').read_text(encoding='utf-8'))
assert protocol['predictor']['seed']==0 and protocol['predictor']['num_diffusion_samples']==1
assert protocol['predictor']['templates'] is False and protocol['predictor']['target_msa_only']
assert sha(root/'inputs/TNFR1_competition_M1_T211.a3m')==protocol['target']['msa_sha256']
if args.archives:
    archive_manifest=json.loads((root/'provenance/archive_manifest.json').read_text(encoding='utf-8'))
    for item in archive_manifest['parts']:
        path=args.archives/item['filename'];assert path.is_file(),path
        assert path.stat().st_size==item['bytes'] and sha(path)==item['sha256'],path
print(f'PASS: {len(rows)} scored candidates; {len(hits)} screening hits; 8 parent trajectories; {verified_files} snapshot files verified.')
