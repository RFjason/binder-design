"""Harvest traceable sequences and build frozen scoring manifests on a CPU node."""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import re
import traceback
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(60, repeat=True)
print('COLLECT_IMPORT', flush=True)

import numpy as np
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
FRAMEWORK = ['QVQLVESGGGLVQPGGSLRLSC', 'WFRQAPGQGLEAVA',
             'YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC', 'WGQGTLVTVS']
assert sum(map(len, FRAMEWORK)) == 83
VHH_PATTERN = re.compile('^' + '([ACDEFGHIKLMNPQRSTVWY]{5,30})'.join(map(re.escape, FRAMEWORK)) + '$')
parser = argparse.ArgumentParser()
parser.add_argument('method', choices=['rfd4','rfantibody','rfantibody_crd3','bc2','bc2_v2'])
parser.add_argument('--batch-id')
parser.add_argument('--campaign', choices=['bulk1','bc2_recovery'],default='bulk1')
args = parser.parse_args()
batch = None
if args.batch_id:
    campaign = json.loads((ROOT/'inputs'/args.campaign/'campaign.json').read_text())
    batch = next(b for b in campaign['batches'] if b['batch_id'] == args.batch_id)
    assert batch['method'] == args.method
label = args.campaign+'_'+args.batch_id if batch else args.method
canonical = json.loads((ROOT / 'P19438.json').read_text())['sequence']['value']
if batch:
    source = Path(batch['output_dir'])
    if args.method == 'rfd4':
        files = sorted((source/'inference_samples').rglob('*_generated.cif'))
    elif args.method == 'rfantibody':
        files = sorted((source/'design_pdbs').glob('*.pdb'))
    else:
        files = sorted(p for p in source.rglob('*') if p.suffix in ('.cif','.pdb') and 'Complexes' in p.parts)
        if args.campaign == 'bc2_recovery':
            files.extend(sorted((source/'1_Trajectories').rglob('trajectory.cif')))
elif args.method == 'rfd4':
    source = ROOT / 'generation/rfd4_canary'
    files = sorted((source / 'inference_samples').rglob('*_generated.cif'))
elif args.method == 'rfantibody':
    source = ROOT / 'generation/rfantibody_canary/design_pdbs'
    files = sorted(source.rglob('*.pdb'))
elif args.method == 'rfantibody_crd3':
    source = ROOT / 'generation/rfantibody_crd3_probe/design_pdbs'
    files = sorted(source.rglob('*.pdb'))
elif args.method == 'bc2':
    source = ROOT / 'generation/bc2_canary_crd2'
    files = sorted([p for p in source.rglob('*') if p.suffix in ('.cif','.pdb')])
else:
    files = sorted(p for folder in (ROOT/'generation').glob('bc2_canary_v2_*')
                   for p in folder.rglob('*') if p.suffix in ('.cif','.pdb') and 'Complexes' in p.parts)
rows, rejects, duplicates, seen = [], [], [], set()
for path in files:
    owned_claim = None
    try:
        print('COLLECT_READ', str(path), flush=True)
        atoms = (get_structure(CIFFile.read(path), model=1)
                 if path.suffix == '.cif' else PDBFile.read(path).get_structure(model=1))
        atoms = atoms[np.isin(atoms.res_name, list(AA)) & (atoms.element != 'H')]
        assert np.isfinite(atoms.coord).all(), 'Non-finite structural coordinates'
        chain_sequences = {str(cid): ''.join(AA[x] for x in atoms.res_name[(atoms.chain_id==cid)&(atoms.atom_name=='CA')])
                           for cid in np.unique(atoms.chain_id)}
        target_chains = [cid for cid, seq in chain_sequences.items() if seq == canonical[29:211]]
        assert len(target_chains) == 1, ('canonical full mature target missing or altered', chain_sequences)
        target_id = target_chains[0]
        binders = [(cid,seq) for cid,seq in chain_sequences.items() if cid != target_id]
        assert len(binders) == 1, chain_sequences
        binder_id, sequence = binders[0]
        assert 100 <= len(sequence) <= 150
        frame_match = VHH_PATTERN.fullmatch(sequence)
        assert frame_match, ('83-residue fixed framework changed or incomplete', sequence)
        assert set(sequence) <= set('ACDEFGHIKLMNPQRSTVWY')
        if sequence in seen:
            continue
        binder = atoms[atoms.chain_id == binder_id]
        sulfur = binder[(binder.res_name == 'CYS') & (binder.atom_name == 'SG')]
        source_scope = 'backbone_only' if set(binder.atom_name) <= {'N','CA','C','O','CB'} else 'all_or_partial_sidechains'
        core_distance, nearest_sg, cysteine_geometry_ok = None, [], None
        if len(sulfur) == sequence.count('C') and len(sulfur) >= 2:
            sd = np.linalg.norm(sulfur.coord[:,None,:] - sulfur.coord[None,:,:], axis=-1)
            np.fill_diagonal(sd, np.inf)
            ca_resids = binder.res_id[binder.atom_name == 'CA']
            core_resids = [ca_resids[len(FRAMEWORK[0])-1], ca_resids[frame_match.start(3)-1]]
            core_indices = [list(sulfur.res_id).index(resid) for resid in core_resids]
            core_distance = float(sd[core_indices[0],core_indices[1]])
            def can_pair(remaining):
                if not remaining: return True
                first, rest = remaining[0], remaining[1:]
                return any(sd[first,partner] <= 3.0 and can_pair([x for x in rest if x != partner]) for partner in rest)
            extras = [index for index in range(len(sulfur)) if index not in core_indices]
            cysteine_geometry_ok = bool(len(sulfur)%2 == 0 and core_distance <= 3.0 and can_pair(extras))
            nearest_sg = sd.min(axis=1).tolist()
        # Keep geometry failures visible as unranked diagnostic canary rows; never call them passing designs.
        if batch:
            assert sequence.count('C') % 2 == 0, 'Odd cysteine count'
            assert cysteine_geometry_ok is not False, 'Source disulfide geometry failed'
        seen.add(sequence)
        sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
        design_id = 'TNFR1_' + label + '_' + sequence_hash[:12]
        if batch:
            registry = ROOT/'candidates/registries/bulk1'
            registry.mkdir(parents=True, exist_ok=True)
            claim = registry/(sequence_hash+'.json')
            try:
                with claim.open('x') as handle:
                    json.dump({'design_id':design_id,'sequence_sha256':sequence_hash,
                               'source_file':str(path),'batch_id':args.batch_id},handle)
                owned_claim = claim
            except FileExistsError:
                try:
                    same_design = json.loads(claim.read_text()).get('design_id') == design_id
                except ValueError:
                    same_design = False
                if not same_design:
                    duplicates.append({'source_file':str(path),'sequence_sha256':sequence_hash,'registry':str(claim)})
                    continue
        output = ROOT / 'candidates' / design_id
        output.mkdir(parents=True, exist_ok=True)
        reference = atoms.copy()
        source_chain_ids = reference.chain_id.copy()
        reference.chain_id[source_chain_ids == target_id] = 'T'
        reference.chain_id[source_chain_ids == binder_id] = 'H'
        # Generator B factors may be absent/NaN. They are not confidence scores.
        # DockQ uses coordinates; preserve the raw source, use neutral metadata in this PDB export.
        reference.set_annotation('b_factor', np.zeros(len(reference), dtype=float))
        pdb = PDBFile(); pdb.set_structure(reference); pdb.write(output / 'designed_complex.pdb')
        (output / 'sequence.fasta').write_text('>'+design_id+'\n'+sequence+'\n')
        row = {'design_id': design_id, 'method': args.method, 'sequence': sequence,
               'structure_method': 'rfantibody' if args.method.startswith('rfantibody') else 'bindcraft2' if args.method.startswith('bc2') else 'rfd4',
               'run_variant': args.method,
               'length': len(sequence), 'source_file': str(path), 'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'canonical_target_verified': True, 'frame_endpoints_checked': True,
               'framework_identity_verified': True, 'fixed_framework_residues': 83,
               'cdr_sequences': list(frame_match.groups()),
               'cdr_lengths': list(map(len, frame_match.groups())),
               'cysteine_geometry_ok': cysteine_geometry_ok, 'cysteine_nearest_SG_A': nearest_sg,
               'cysteine_count': sequence.count('C'), 'source_SG_atom_count': len(sulfur),
               'source_atom_scope': source_scope,
               'source_disulfide_status': 'ASSESSED' if cysteine_geometry_ok is not None else 'NOT_ASSESSABLE_WITH_MISSING_SIDECHAINS',
               'core_cysteine_SG_distance_A': core_distance,
               'reference_pdb': str(output / 'designed_complex.pdb'),
               'reference_pdb_b_factor': 'zero placeholder; never used as predicted confidence',
               'status': 'CANARY_DIAGNOSTIC_NOT_FINAL_SUBMISSION',
               'novelty_status': 'NOT_RUN',
               'vhh_numbering_status': 'EXACT_SOURCE_FRAMEWORK_AND_THREE_CDR_SEGMENTS_VERIFIED'}
        if batch:
            stem = path.stem
            parent = stem.split('_dldesign_')[0] if args.method == 'rfantibody' else stem.split('_candidate')[0] if args.method == 'bc2_v2' else stem
            row.update(campaign=args.campaign,batch_id=args.batch_id,stratum=batch['stratum'],
                       parent_id=args.batch_id+':'+('trajectory_1' if args.method=='bc2_v2' else parent),generation_seed=batch['seed'],
                       status='BULK_CANDIDATE_PENDING_SCORE',run_variant=args.campaign+'_'+args.method,
                       source_generation_stage='gradient_trajectory' if path.name=='trajectory.cif' else 'native_generated_or_refolded')
        (output / 'provenance.json').write_text(json.dumps(row, indent=2))
        rows.append(row)
    except Exception as exc:
        if owned_claim is not None and owned_claim.exists():
            try:
                if not any(r['design_id']==design_id for r in rows):
                    owned_claim.unlink()
            except (OSError, ValueError):
                pass
        rejects.append({'source': str(path), 'reason': repr(exc), 'traceback': traceback.format_exc()})
cases = []
for row in rows:
    cases.extend([
        {'id': row['design_id'] + '_complex', 'candidate_id': row['design_id'],
         'chains': [{'id': 'A', 'sequence': canonical[:211], 'msa': 'inputs/TNFR1_competition_M1_T211.a3m'},
                    {'id': 'B', 'sequence': row['sequence']}]},
        {'id': row['design_id'] + '_monomer', 'candidate_id': row['design_id'],
         'chains': [{'id': 'B', 'sequence': row['sequence']}]},
    ])
manifest = ROOT / 'inputs' / ('canary_scoring_' + label + '.json')
manifest.write_text(json.dumps(cases, indent=2))
report = {'method': args.method, 'source_files': len(files), 'candidates': rows, 'rejected_files': rejects,
          'duplicates_skipped':duplicates,'batch':batch,
          'scoring_manifest': str(manifest), 'status': 'COLLECTED' if rows else 'NO_VALID_CANDIDATES'}
(ROOT / 'reports' / ('canary_collection_' + label + '.json')).write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2), flush=True)
faulthandler.cancel_dump_traceback_later()
if not rows and not batch:
    raise SystemExit('No structurally complete candidates from '+args.method)
