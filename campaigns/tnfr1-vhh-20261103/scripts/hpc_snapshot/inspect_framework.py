"""Inspect official VHH conditions inside the installed runtime, on Slurm."""
import json
import os
from pathlib import Path
import numpy as np
import rfproteina.conditions
from atomworks.io import parse
from atomworks.io.utils.atom_array_plus import as_atom_array_plus
from rfproteina.utils.parse_presets import ANNOTATED_CIF

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
SOURCE = Path('/dssg/home/acct-clshzh/jasoncao1280/rfd4_deploy_20261022/src/RFD4-Proteina')
path = SOURCE / 'examples/antibody_design/structures/vhh/pmvk_3CH4__h-NbBCII10_cdr123.cif.zst'
aa = as_atom_array_plus(parse(path, config=ANNOTATED_CIF)['assemblies']['1'][0])
summary = {'annotations': aa.get_annotation_categories(), 'chains': {}}
for cid in np.unique(aa.chain_id):
    chain = aa[aa.chain_id == cid]
    summary['chains'][cid] = {'atoms': len(chain), 'resids': np.unique(chain.res_id).tolist()}
for annotation in aa.get_annotation_categories():
    if 'condition' in annotation or 'expand' in annotation:
        value = aa.get_annotation(annotation)
        summary[annotation] = str(value[:5]) if value.ndim > 1 else list(map(str, np.unique(value)))[:30]
print(json.dumps(summary, indent=2), flush=True)
(ROOT / 'reports/framework_annotation_inspection.json').write_text(json.dumps(summary, indent=2))
from atomworks.io.utils.io_utils import to_cif_file
to_cif_file(aa, ROOT / 'sources/official_vhh_example.cif', save_standard_annotations=True, chain_disambiguation='chain_iid')
