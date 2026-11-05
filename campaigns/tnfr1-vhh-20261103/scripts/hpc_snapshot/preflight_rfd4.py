"""Run the installed model's input featurization/validation without loading weights."""
import json
import os
from pathlib import Path
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(60, repeat=True)
print('PREFLIGHT_IMPORT_TORCH', flush=True)
import torch
print('PREFLIGHT_IMPORT_PIPELINE', flush=True)
from rfproteina.datapipes.pipeline import build_design_validation_pipeline
from rfproteina.inference_utils.datasets import DesignProblemDataset

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
print('PREFLIGHT_BUILD_PIPELINE', flush=True)
pipeline = build_design_validation_pipeline(metrics_tags=[], sampling_transforms=None)
print('PREFLIGHT_BUILD_DATASET', flush=True)
dataset = DesignProblemDataset(dir_path_or_file_paths=str(ROOT/'inputs/rfd4'), transform=pipeline,
                               num_replicates=1, prevalidate=True)
rows = []
for index in range(len(dataset)):
    print('PREFLIGHT_ITEM', index, flush=True)
    data = dataset[index]
    assert data is not None
    rows.append({'index': index, 'keys': list(data),
                 'tensor_shapes': {key: list(value.shape) for key,value in data.items() if isinstance(value,torch.Tensor)}})
result = {'ok': True, 'scope': 'official CPU input validation and featurization; no weights or inference',
          'job_id': os.environ['SLURM_JOB_ID'], 'inputs': rows}
(ROOT/'reports/rfd4_input_preflight.json').write_text(json.dumps(result, indent=2))
print('RFD4_INPUT_PREFLIGHT_OK', json.dumps(result), flush=True)
faulthandler.cancel_dump_traceback_later()
