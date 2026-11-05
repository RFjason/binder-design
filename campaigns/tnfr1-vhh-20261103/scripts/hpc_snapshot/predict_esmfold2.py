"""Frozen single-seed ESMFold2 inference, on a Slurm compute node only."""
import argparse
import hashlib
import importlib.metadata as md
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import torch
from esm.models.esmfold2 import ESMFold2InputBuilder, EsmFold2Model, ProteinInput, StructurePredictionInput
from esm.utils.msa import MSA

if not os.environ.get('SLURM_JOB_ID'):
    raise SystemExit('Slurm compute nodes only')
ROOT = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
LEGACY = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006')
KWARGS = dict(num_loops=20, num_sampling_steps=200, num_diffusion_samples=1,
              seed=0, msa_max_depth=1024, msa_column_mask_rate=0.1, lm_dropout=0.3)
protocol = json.loads((ROOT / 'scoring_protocol_v1.json').read_text())
assert all(protocol['predictor'][key] == value for key,value in KWARGS.items())
canonical = json.loads((ROOT / 'P19438.json').read_text())['sequence']['value']
assert hashlib.sha256(canonical[:211].encode()).hexdigest() == protocol['target']['competition_sequence_sha256']
assert hashlib.sha256((ROOT / protocol['target']['msa']).read_bytes()).hexdigest() == protocol['target']['msa_sha256']
assert hashlib.sha256((LEGACY / 'src/IPSAE/ipsae.py').read_bytes()).hexdigest() == protocol['score']['reference_script_sha256']
parser = argparse.ArgumentParser()
parser.add_argument('manifest')
args = parser.parse_args()
cases = json.loads(Path(args.manifest).read_text())
model = EsmFold2Model.from_pretrained(str(LEGACY / 'weights/ESMFold2'), device='cuda').eval()
builder = ESMFold2InputBuilder()
failures = []
for item in cases:
    out = ROOT / 'predictions' / item['id']
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'summary.json').exists() and json.loads((out / 'summary.json').read_text()).get('ok'):
        continue
    try:
        proteins, chain_ids, sequences, msa_records = [], [], [], []
        for entry in item['chains']:
            sequence = entry['sequence']
            assert set(sequence) <= set('ACDEFGHIKLMNPQRSTVWY')
            msa = MSA.from_a3m(ROOT / entry['msa'], remove_insertions=True) if entry.get('msa') else None
            if msa is not None:
                assert msa.sequences[0] == sequence and msa.depth > 1
            proteins.append(ProteinInput(id=entry['id'], sequence=sequence, msa=msa))
            chain_ids.append(entry['id'])
            sequences.append(sequence)
            msa_records.append({'chain': entry['id'], 'msa': entry.get('msa'), 'rows': msa.depth if msa else 0,
                                'sha256': hashlib.sha256((ROOT / entry['msa']).read_bytes()).hexdigest() if msa else None})
        case = StructurePredictionInput(sequences=proteins)
        features, chains = builder.prepare_input(case, seed=0)
        indices, token_chain_ids, residue_numbers = [], [], []
        for chain in chains:
            for token in chain.tokens:
                indices.append(token.token_index)
                token_chain_ids.append(chain.chain_id)
                residue_numbers.append(token.residue_index)
        count = sum(map(len, sequences))
        assert len(indices) == count and sorted(indices) == list(range(count))
        for cid, seq in zip(chain_ids, sequences):
            assert token_chain_ids.count(cid) == len(seq)
        del features
        (out / 'input_and_token_mapping.json').write_text(json.dumps({
            'case': item, 'msa_records': msa_records, 'kwargs': KWARGS, 'templates': False,
            'initial_guess': False, 'cif_order_to_pae_token_index': indices,
            'cif_token_chain_ids': token_chain_ids, 'cif_token_residue_numbers': residue_numbers}, indent=2))
        torch.cuda.reset_peak_memory_stats()
        start = time.monotonic()
        print('START_CASE', item['id'], flush=True)
        with torch.inference_mode():
            result = builder.fold(model, case, complex_id=item['id'], **KWARGS)
        torch.cuda.synchronize()
        elapsed = time.monotonic() - start
        raw_pae = torch.as_tensor(result.pae).float().cpu().numpy()
        raw_plddt = torch.as_tensor(result.plddt).float().cpu().numpy()
        assert raw_pae.shape == (count, count) and raw_plddt.shape == (count,)
        assert np.isfinite(raw_pae).all() and np.isfinite(raw_plddt).all()
        pae = raw_pae[np.ix_(indices, indices)]
        plddt = raw_plddt[indices]
        cif = result.complex.to_mmcif()
        (out / 'model.cif').write_text(cif)
        np.savez_compressed(out / 'pae.npz', pae=raw_pae, cif_order_pae=pae, plddt=raw_plddt,
                            cif_order_to_pae_token_index=np.array(indices), cif_token_chain_ids=np.array(token_chain_ids))
        fields, atom_plddt, ca_chains = [], [], []
        for line in cif.splitlines():
            if line.startswith('_atom_site.'):
                fields.append(line.split('.', 1)[1].strip())
            elif line.startswith(('ATOM ', 'HETATM ')):
                values = line.split()
                atom_plddt.append(float(values[fields.index('B_iso_or_equiv')]))
                if values[fields.index('label_atom_id')].strip('"') == 'CA':
                    ca_chains.append(values[fields.index('label_asym_id')])
        assert ca_chains == token_chain_ids, 'CIF/PAE mapping mismatch'
        scale = 100 if plddt.max() <= 1.1 else 1
        summary = {'ok': True, 'id': item['id'], 'job_id': os.environ['SLURM_JOB_ID'],
                   'protocol_id': protocol['protocol_id'],
                   'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'hostname': os.uname().nodename, 'esm': md.version('esm'),
                   'kwargs': KWARGS, 'templates': False, 'initial_guess': False,
                   'msa_records': msa_records, 'fold_seconds': elapsed,
                   'peak_gpu_memory_bytes': torch.cuda.max_memory_allocated(),
                   'pae_cif_mapping_validated': True,
                   'chain_plddt': {cid: float(plddt[np.array(token_chain_ids) == cid].mean() * scale) for cid in chain_ids}}
        if len(chain_ids) == 2:
            if max(atom_plddt) <= 1.1:
                atom_plddt = [v * 100 for v in atom_plddt]
            (out / 'af3like_pae.json').write_text(json.dumps({'pae': pae.tolist(), 'atom_plddts': atom_plddt,
                 'token_chain_ids': token_chain_ids, 'token_res_ids': residue_numbers}))
            script = LEGACY / 'src/IPSAE/ipsae.py'
            process = subprocess.run([sys.executable, str(script), str(out / 'af3like_pae.json'),
                                      str(out / 'model.cif'), '10', '10'], capture_output=True, text=True)
            (out / 'ipsae_stdout.txt').write_text(process.stdout)
            (out / 'ipsae_stderr.txt').write_text(process.stderr)
            assert process.returncode == 0, process.stderr
            directions = {}
            for line in (out / 'model_10_10.txt').read_text().splitlines():
                values = line.split()
                if len(values) > 6 and values[4] == 'asym':
                    directions[values[0] + '_to_' + values[1]] = float(values[5])
            assert len(directions) == 2
            summary.update(ipsae_directions=directions, ipsae_min=min(directions.values()),
                           ipsae_max=max(directions.values()), ipTM=float(result.iptm),
                           ipsae_script_sha256=hashlib.sha256(script.read_bytes()).hexdigest())
        from biotite.structure.io.pdbx import CIFFile, get_structure
        from biotite.structure.io.pdb import PDBFile
        atoms = get_structure(CIFFile.read(io.StringIO(cif)), model=1, extra_fields=['b_factor'])
        pdb = PDBFile()
        pdb.set_structure(atoms)
        pdb.write(out / 'model.pdb')
        (out / 'summary.json').write_text(json.dumps(summary, indent=2))
        print('CASE_OK', json.dumps(summary), flush=True)
        del result
    except Exception as exc:
        import traceback
        failure = {'ok': False, 'id': item['id'], 'error': repr(exc), 'traceback': traceback.format_exc()}
        (out / ('failure_' + os.environ['SLURM_JOB_ID'] + '.json')).write_text(json.dumps(failure, indent=2))
        print('CASE_FAIL', json.dumps(failure), flush=True)
        failures.append(item['id'])
if failures:
    raise SystemExit('Failed cases: ' + ','.join(failures))
