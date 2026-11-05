"""Close the canary loop with frozen-score checks and pose consistency."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import re

import numpy as np
from biotite.structure.io.pdb import PDBFile

assert os.environ.get('SLURM_JOB_ID'), 'Slurm only'
ROOT = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103')
LEGACY = Path('/dssg/home/acct-clshzh/jasoncao1280/benchmarks/anthropic_pdl1_pipeline_20261006')
parser = argparse.ArgumentParser(); parser.add_argument('method'); parser.add_argument('--batch-id')
parser.add_argument('--campaign',choices=['bulk1','bc2_recovery'],default='bulk1');args = parser.parse_args()
label = args.campaign+'_'+args.batch_id if args.batch_id else args.method
collection = json.loads((ROOT / 'reports' / ('canary_collection_'+label+'.json')).read_text())
manifest = json.loads((ROOT / 'inputs' / ('canary_scoring_'+label+'.json')).read_text())
requested_ids = {case['candidate_id'] for case in manifest}
rows, failures = [], []
for candidate in collection['candidates']:
    cid = candidate['design_id']
    if cid not in requested_ids:
        continue
    try:
        complex_folder = ROOT / 'predictions' / (cid+'_complex')
        monomer_folder = ROOT / 'predictions' / (cid+'_monomer')
        score = json.loads((complex_folder/'summary.json').read_text())
        monomer = json.loads((monomer_folder/'summary.json').read_text())
        assert score['ok'] and monomer['ok'] and score['pae_cif_mapping_validated']
        assert score['kwargs']['seed'] == 0 and score['kwargs']['num_diffusion_samples'] == 1
        assert score['templates'] is False and score['initial_guess'] is False
        assert score['msa_records'][0]['rows'] > 1 and score['msa_records'][1]['rows'] == 0
        predicted = PDBFile.read(complex_folder/'model.pdb').get_structure(model=1)
        original_ids = predicted.chain_id.copy()
        # Build the mature target subset by sequence order, rather than assuming a token-index convention.
        target_ca = predicted[(original_ids == 'A') & (predicted.atom_name == 'CA')]
        binder_atoms=predicted[original_ids=='B']
        signal_atoms=predicted[(original_ids=='A') & np.isin(predicted.res_id,target_ca.res_id[:29])]
        signal_distances=np.linalg.norm(binder_atoms.coord[:,None,:]-signal_atoms.coord[None,:,:],axis=-1)
        signal_contact_resids=np.unique(binder_atoms.res_id[signal_distances.min(axis=1)<=4.5])
        mature_resids = target_ca.res_id[29:]
        mature = predicted[((original_ids=='A') & np.isin(predicted.res_id,mature_resids)) | (original_ids=='B')]
        mature.chain_id[mature.chain_id=='A'] = 'T'
        mature.chain_id[mature.chain_id=='B'] = 'H'
        pdb = PDBFile(); pdb.set_structure(mature); pdb.write(complex_folder/'mature_pose.pdb')
        dockq_json = complex_folder/'self_consistency_dockq.json'
        call = subprocess.run([str(LEGACY/'envs/scoring/bin/DockQ'), str(complex_folder/'mature_pose.pdb'),
                               candidate['reference_pdb'], '--mapping', 'TH:TH', '--n_cpu', '1',
                               '--json', str(dockq_json)], capture_output=True, text=True, timeout=180)
        (complex_folder/'dockq_stdout.txt').write_text(call.stdout)
        (complex_folder/'dockq_stderr.txt').write_text(call.stderr)
        pose_status = 'ASSESSED'
        source_interface_diagnostic = None
        if call.returncode != 0 and 'Could not find interfaces in the native model' in call.stderr:
            # Preserve completed primary predictions when the generated reference has
            # no DockQ-assessable interface. A missing pose metric cannot pass selection.
            ref = PDBFile.read(candidate['reference_pdb']).get_structure(model=1)
            rt = ref[ref.chain_id=='T']; rh = ref[ref.chain_id=='H']
            rd = np.linalg.norm(rt.coord[:,None,:]-rh.coord[None,:,:],axis=-1)
            source_interface_diagnostic = {'minimum_atom_distance_A':float(rd.min()),
                                           'atom_pairs_within_5A':int((rd<5.0).sum())}
            pose_dockq = None
            pose_status = 'NO_DOCKQ_ASSESSABLE_INTERFACE_IN_DESIGN_REFERENCE'
        else:
            assert call.returncode == 0 and dockq_json.exists(), call.stderr
            dockq = json.loads(dockq_json.read_text())
            assert len(dockq['best_result']) == 1
            pose_dockq = float(dockq['GlobalDockQ'])
        single_plddt = float(monomer['chain_plddt']['B'])
        mono_atoms = PDBFile.read(monomer_folder/'model.pdb').get_structure(model=1)
        mono_ca = mono_atoms[mono_atoms.atom_name=='CA']
        mono_sg = mono_atoms[(mono_atoms.res_name=='CYS') & (mono_atoms.atom_name=='SG')]
        sg_dist = np.linalg.norm(mono_sg.coord[:,None,:]-mono_sg.coord[None,:,:],axis=-1)
        np.fill_diagonal(sg_dist,np.inf)
        source_frame = ['QVQLVESGGGLVQPGGSLRLSC','WFRQAPGQGLEAVA','YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC','WGQGTLVTVS']
        frame = re.fullmatch('([ACDEFGHIKLMNPQRSTVWY]{5,30})'.join(map(re.escape,source_frame)),candidate['sequence'])
        assert frame and len(mono_ca)==candidate['length']
        assert len(mono_sg)==candidate['sequence'].count('C')
        core_resids=[mono_ca.res_id[len(source_frame[0])-1],mono_ca.res_id[frame.start(3)-1]]
        core_indices=[list(mono_sg.res_id).index(r) for r in core_resids]
        core_distance=float(sg_dist[core_indices[0],core_indices[1]])
        def pairable(indices):
            if not indices: return True
            first,rest=indices[0],indices[1:]
            return any(sg_dist[first,other]<=3.0 and pairable([x for x in rest if x!=other]) for other in rest)
        extras=[i for i in range(len(mono_sg)) if i not in core_indices]
        monomer_disulfides_ok=bool(len(mono_sg)%2==0 and core_distance<=3.0 and pairable(extras))
        ectodomain=predicted[(original_ids=='A') & np.isin(predicted.res_id,mature_resids)]
        contact_distances=np.linalg.norm(binder_atoms.coord[:,None,:]-ectodomain.coord[None,:,:],axis=-1)
        interface_resids=np.unique(binder_atoms.res_id[contact_distances.min(axis=1)<=4.5])
        bca=binder_atoms[binder_atoms.atom_name=='CA']
        interface_mask=np.isin(bca.res_id,interface_resids)
        full_cdr_mask=np.zeros(candidate['length'],dtype=bool)
        for group in (1,2,3): full_cdr_mask[frame.start(group):frame.end(group)]=True
        cdr_fraction=float(full_cdr_mask[interface_mask].mean()) if interface_mask.any() else None
        glycan_clearance={}
        for pos in (54,145,151):
            resid=target_ca.res_id[pos-1]
            anchor=predicted[(original_ids=='A') & (predicted.res_id==resid) & (predicted.atom_name=='ND2')]
            glycan_clearance[str(pos)]=float(np.linalg.norm(binder_atoms.coord-anchor.coord[0],axis=1).min()) if len(anchor) else None
        row = {**candidate, 'ipsae_min': score['ipsae_min'], 'ipsae_directions': score['ipsae_directions'],
               'ipTM': score['ipTM'], 'complex_binder_plddt': score['chain_plddt']['B'],
               'monomer_plddt': single_plddt, 'pose_dockq': pose_dockq,
               'pose_dockq_status':pose_status,'source_interface_diagnostic':source_interface_diagnostic,
               'monomer_disulfide_geometry_ok':monomer_disulfides_ok,
               'monomer_core_SG_distance_A':core_distance,
               'predicted_full_CDR_interface_fraction_4p5A':cdr_fraction,
               'predicted_signal_peptide_contact_binder_residues_4p5A':len(signal_contact_resids),
               'mature_construct_validation_required':bool(len(signal_contact_resids)),
               'predicted_glycan_anchor_clearance_A':glycan_clearance,
               'scoring_protocol_verified': True, 'pipeline_canary_pass': True,
               'computational_quality_pass': bool(single_plddt >= 70 and pose_dockq is not None and pose_dockq >= 0.23
                        and monomer_disulfides_ok and score['ipsae_min'] >= 0.6),
               'final_rank_eligible': False, 'pending_final_gate': 'novelty/developability/selection audit'}
        rows.append(row)
    except Exception as exc:
        failures.append({'design_id': cid, 'error': repr(exc)})
result = {'method': args.method, 'batch_id':args.batch_id,
          'pipeline_pass': (bool(rows) or (bool(args.batch_id) and not requested_ids)) and not failures,
          'total_collected': len(collection['candidates']), 'requested_scoring_count': len(requested_ids),
          'computational_quality_pass_count': sum(r['computational_quality_pass'] for r in rows),
          'scored_candidates': rows, 'failures': failures,
          'scope': 'canary loop validation, not experimental binding evidence'}
(ROOT/'reports'/('canary_result_'+label+'.json')).write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2), flush=True)
if not result['pipeline_pass']:
    raise SystemExit('Canary scoring/pose adapter incomplete')
