"""Prepare complete-ectodomain, constrained VHH inputs inside Slurm."""
import hashlib
import itertools
import json
import os
from pathlib import Path

import numpy as np
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile
import rfproteina.conditions
from atomworks.io import parse
from atomworks.io.utils.atom_array_plus import as_atom_array_plus, concatenate_atom_array_plus
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.ml.conditions import C_SEQ, C_CRD
from rfproteina.conditions import C_HOT, C_CRDGRP, C_CATH, C_INTCRD
from rfproteina.utils.parse_presets import ANNOTATED_CIF, DISTILLATION

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
SOURCE = Path('/dssg/home/acct-clshzh/jasoncao1280/rfd4_deploy_20261022/src/RFD4-Proteina')
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
canonical = json.loads((ROOT / 'P19438.json').read_text())['sequence']['value']
prediction = ROOT / 'inputs/TNFR1_AFDB_mature_L30_T211.pdb'
PDB_PREDICTION = DISTILLATION.replace(load_standard_annotations=False)
target = as_atom_array_plus(parse(prediction, config=PDB_PREDICTION)['assemblies']['1'][0])
target = target[(target.chain_id == 'A') & np.isin(target.res_name, list(AA))]
assert np.isfinite(target.coord).all()
for residue in struc.residue_iter(target):
    assert set(['N','CA','C','O']) <= set(residue.atom_name)
assert ''.join(AA[x] for x in target.res_name[target.atom_name == 'CA']) == canonical[29:211]
# AFDB input already carries canonical residue numbering.
assert target.res_id[target.atom_name == 'CA'].tolist() == list(range(30, 212))
C_SEQ.set_annotation_from_ground_truth(target, np.ones(len(target), dtype=bool))
backbone = np.isin(target.atom_name, ['N', 'CA', 'C', 'O'])
C_CRD.set_annotation_from_ground_truth(target, backbone)
C_HOT.set_annotation(target, np.zeros(len(target), dtype=bool))

input_dir = ROOT / 'inputs'
input_dir.mkdir(exist_ok=True)
to_cif_file(target, input_dir / 'TNFR1_mature_canonical.cif', save_standard_annotations=True, chain_disambiguation='chain_iid')
pdb = PDBFile(); pdb.set_structure(target); pdb.write(input_dir / 'TNFR1_mature_canonical.pdb')
sasa = struc.sasa(target, point_number=500)
glycan_positions = [54, 145, 151]
anchors = target.coord[np.isin(target.res_id, glycan_positions) & (target.atom_name == 'ND2')]
assert len(anchors) == 3
ca = target[target.atom_name == 'CA']
experimental = PDBFile.read(ROOT / '9S3I.pdb').get_structure(model=1)
experimental = experimental[(experimental.chain_id == 'A') & (experimental.atom_name == 'CA')]
experimental.res_id += 29
domain_rmsds = {}
for start, end in [(42,81), (84,125), (127,166)]:
    observed = experimental[(experimental.res_id >= start) & (experimental.res_id <= end)]
    pred = ca[np.isin(ca.res_id, observed.res_id)]
    assert np.array_equal(pred.res_id, observed.res_id)
    fitted, _ = struc.superimpose(observed.coord, pred.coord)
    domain_rmsds[f'{start}-{end}'] = float(struc.rmsd(observed.coord, fitted))
assert max(domain_rmsds.values()) < 3.5, domain_rmsds

patches = {'crd2': [96, 98, 100, 104, 106, 108, 109],
           'crd3': [129, 132, 133, 134, 135, 136, 138, 140, 141, 142]}
report = {'job_id': os.environ['SLURM_JOB_ID'], 'target_length': 182,
          'design_context_source': 'AFDB AF-P19438-F1-model_v6, canonical mature L30-T211',
          'design_context_sha256': hashlib.sha256(prediction.read_bytes()).hexdigest(),
          'target_sequence_sha256': hashlib.sha256(canonical[29:211].encode()).hexdigest(),
          'experimental_domain_CA_RMSD_A': domain_rmsds,
          'potential_glycan_anchors': glycan_positions, 'patches': {}, 'inputs': []}
template = as_atom_array_plus(parse(SOURCE / 'examples/antibody_design/structures/vhh/pmvk_3CH4__h-NbBCII10_cdr123.cif.zst',
                                   config=ANNOTATED_CIF)['assemblies']['1'][0])
framework = template[template.chain_id == 'H']
assert np.unique(framework.res_id).size == 86
for patch_name, positions in patches.items():
    candidate_rows = []
    for position in positions:
        mask = (target.res_id == position) & ~np.isin(target.atom_name, ['N','CA','C','O'])
        candidate_indices = np.flatnonzero(mask & (sasa >= 1.0))
        if not len(candidate_indices):
            continue
        distances = np.sqrt(((target.coord[candidate_indices,None,:] - anchors[None,:,:])**2).sum(axis=-1)).min(axis=1)
        eligible = candidate_indices[distances >= 8.0]
        if not len(eligible):
            continue
        index = int(eligible[np.argmax(sasa[eligible])])
        candidate_rows.append({'position': position, 'aa': canonical[position-1], 'atom': str(target.atom_name[index]),
                               'sasa_A2': float(sasa[index]), 'glycan_anchor_distance_A': float(np.linalg.norm(target.coord[index]-anchors, axis=1).min()),
                               'index': index})
    candidate_rows.sort(key=lambda x: x['sasa_A2'], reverse=True)
    # Pick a compact patch jointly; the most exposed residue alone may be isolated.
    selected = []
    for count in (3, 2):
        combinations = [group for group in itertools.combinations(candidate_rows, count)
                        if all(np.linalg.norm(target.coord[a['index']]-target.coord[b['index']]) <= 15.0
                               for a,b in itertools.combinations(group,2))]
        if combinations:
            selected = list(max(combinations, key=lambda group: sum(x['sasa_A2'] for x in group)))
            break
    assert len(selected) >= 2, (patch_name, candidate_rows)
    hotmask = np.zeros(len(target), dtype=bool)
    hotmask[[x['index'] for x in selected]] = True
    report['patches'][patch_name] = {'selected_hotspot_atoms': selected, 'all_considered': candidate_rows,
         'interpretation': 'Solvent-exposed structural candidate patch, not experimentally validated energetic hotspots'}
    for length_name, ranges in {'compact': [(11,13),(8,10),(10,14)],
                                 'long': [(13,15),(9,11),(16,22)]}.items():
        receptor = target.copy()
        C_HOT.set_annotation(receptor, hotmask)
        binder = framework.copy()
        placeholder_ids = np.unique(binder.res_id[binder.get_annotation('annotation_expsegmin_1_residue') >= 0])
        assert len(placeholder_ids) == 3
        for resid, (lo, hi) in zip(placeholder_ids, ranges):
            mask = binder.res_id == resid
            binder.set_annotation('annotation_expsegmin_1_residue', np.where(mask, lo, binder.get_annotation('annotation_expsegmin_1_residue')))
            binder.set_annotation('annotation_expsegmax_1_residue', np.where(mask, hi, binder.get_annotation('annotation_expsegmax_1_residue')))
            # Keep the source placeholder CCD identifier/registry; sampling bounds live in S_SEGMIN/MAX.
        # Strip irrelevant hotspots left on the source binder, preserve its fixed framework/CATH conditions.
        C_HOT.set_annotation(binder, np.zeros(len(binder), dtype=bool))
        combined = concatenate_atom_array_plus([receptor, binder], fill_missing_standard_annotations=True,
                                               on_annotation_mismatch_policy='warn')
        output = input_dir / 'rfd4' / f'tnfr1_{patch_name}_{length_name}.cif'
        output.parent.mkdir(exist_ok=True)
        to_cif_file(combined, output, save_standard_annotations=True, chain_disambiguation='chain_iid')
        reread = as_atom_array_plus(parse(output, config=ANNOTATED_CIF)['assemblies']['1'][0])
        assert set(reread.chain_id) == {'A','H'}
        assert len(np.unique(reread.res_id[reread.chain_id == 'A'])) == 182
        assert int(C_HOT.annotation(reread).sum()) == len(selected)
        saved_binder = reread[reread.chain_id == 'H']
        saved_placeholders = np.unique(saved_binder.res_id[saved_binder.get_annotation('annotation_expsegmin_1_residue') >= 0])
        assert len(saved_placeholders) == 3, 'Variable CDR annotations lost in concatenation/write'
        for resid, (lo,hi) in zip(saved_placeholders,ranges):
            mask = saved_binder.res_id == resid
            assert np.all(saved_binder.get_annotation('annotation_expsegmin_1_residue')[mask] == lo)
            assert np.all(saved_binder.get_annotation('annotation_expsegmax_1_residue')[mask] == hi)
        assert any(str(x).startswith('2.60.40') for x in C_CATH.annotation(saved_binder))
        assert np.isfinite(C_INTCRD.annotation(saved_binder)).all(axis=-1).sum() >= 300
        report['inputs'].append({'path': str(output), 'sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
                                  'framework_residues': 83, 'CDR_ranges': ranges,
                                  'minimum_vhh_length': 83+sum(x[0] for x in ranges),
                                  'maximum_vhh_length': 83+sum(x[1] for x in ranges),
                                  'cdr_length_CATH_internal_coordinate_annotations_roundtrip_verified': True})
    # Full canonical-numbered receptor shared with RFantibody/BindCraft2.
    (input_dir / f'hotspots_{patch_name}.txt').write_text(','.join('A'+str(x['position']) for x in selected)+'\n')
report['sequence_conditions'] = 'target and framework fixed; only CDR placeholder sequence/structure generated'
# BindCraft2 takes the independently sourced, complete HLT framework, not its bundled historical design.
rfab_root = Path('/dssg/home/acct-clshzh/jasoncao1280/RFantibody')
raw_framework = as_atom_array_plus(parse(rfab_root / 'scripts/examples/example_inputs/h-NbBCII10.pdb',
                                        config=PDB_PREDICTION)['assemblies']['1'][0])
raw_framework = raw_framework[raw_framework.chain_id == 'H']
scaffold_file = input_dir / 'h-NbBCII10_complete.cif'
to_cif_file(raw_framework, scaffold_file, save_standard_annotations=True, chain_disambiguation='chain_iid')
framework_seq = ''.join(AA[x] for x in raw_framework.res_name[raw_framework.atom_name == 'CA'])
assert 100 <= len(framework_seq) <= 150
report['bc2_framework'] = {'sequence': framework_seq, 'source': str(rfab_root / 'scripts/examples/example_inputs/h-NbBCII10.pdb'),
                            'sha256': hashlib.sha256(scaffold_file.read_bytes()).hexdigest()}
for patch_name in patches:
    config = {'core': 'benchmark', 'modality': 'VHH',
              'targets': [{'name': 'TNFR1', 'target_path': str(input_dir / 'TNFR1_mature_canonical.pdb'),
                           'chains': 'A', 'hotspots': ','.join(str(x['position']) for x in report['patches'][patch_name]['selected_hotspot_atoms'])}],
              'binder_scaffold': str(scaffold_file),
              'mutate_positions': 'H26-35,H55-60,H102-117',
              'weights_binder_intra_coldspot': 0.0, 'weights_binder_coldspot': 0.0,
              'min_scaffold_sequence_retained_final': 1.0,
              'number_of_final_designs': 1, 'max_trajectories': 1,
              'auto_multi_gpu': False, 'design_workers': 1,
              'min_plddt_refine': 0.65, 'min_plddt_mutate': 0.65,
              'min_iptm_final': 0.6, 'min_monomer_plddt_final': 0.7,
              'project_folder': str(ROOT / 'generation' / ('bc2_canary_' + patch_name)),
              'save_failed_refolds': True, 'save_failed_trajectories': True,
              'archive_trajectories': False, 'campaign_seed': 301 if patch_name == 'crd2' else 302}
    (input_dir / ('bc2_canary_' + patch_name + '.json')).write_text(json.dumps(config, indent=2))
(ROOT / 'reports/design_input_audit.json').write_text(json.dumps(report, indent=2))
print('DESIGN_INPUTS_OK', json.dumps(report), flush=True)
