"""Validate parent loop masks using the native RFantibody parser on CPU Slurm."""
import os
import json
from pathlib import Path
from rfantibody.util.pose import Pose
from rfantibody.proteinmpnn.sample_features import SampleFeatures
assert os.environ.get('SLURM_JOB_ID')
root=Path(__file__).resolve().parents[1]
for path in sorted((root/'inputs/parents').glob('*/parent.pdb')):
 pose=Pose.from_pdb(str(path))
 print(path,vars(pose).keys(),flush=True)
 assert str(pose.chain[0])=='H' and str(pose.chain[-1])=='T'
 features=SampleFeatures(pose,str(path));features.loop_string2fixed_res('H1,H2,H3')
 assert all(pose.cdr_dict[k] for k in ('H1','H2','H3'))
 assert len(features.fixed_res['T'])==182
 # Native annotation accessor is checked in the preflight before jobs release.
 assert sum(' H1' in line for line in path.read_text().splitlines())>0
 assert sum(' H2' in line for line in path.read_text().splitlines())>0
 assert sum(' H3' in line for line in path.read_text().splitlines())>0
(root/'reports/parent_native_ok.json').write_text(json.dumps({'ok':True,'job_id':os.environ['SLURM_JOB_ID']}))
