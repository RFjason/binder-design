"""Record recovery IDs without overwriting the failed-run history."""
import json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
path=root/'reports/job_ledger.json'
ledger=json.loads(path.read_text())
if 'initial_routes' not in ledger:
    ledger['initial_routes']=ledger['routes']
ledger['routes']={
 'rfd4':{'generation':'64563125','collection':'64563167','prediction':'64563170','audit':'64563174'},
 'rfantibody':{'generation':'64563119','collection':'64563168','prediction':'64563171','audit':'64563175'},
 'bc2':{'generation':'64536750','collection':'64536759','prediction':'64536763','audit':'64536772'},
 'bc2_v2':{'generation_array':'64563169','collection':'64563172','prediction':'64563176','audit':'64563179'}}
ledger['recovery_preflights']={'rfd4_launcher':'64563118','bc2_v2':'64563164'}
ledger['first_generation_failures']={
 '64536748':'Hydra override ordering in shared launcher; failed before generation',
 '64536749':'RFantibody subprocess could not find python; failed before generation'}
ledger['observed_results']={'target_prediction':'COMPLETED', 'bc2_v1_unique_sequences':10,
 'bc2_v1_length':127,'bc2_v1_internal_accepted':0,'bc2_v1_cysteine_count':5,
 'bc2_v1_cysteine_geometry_pass':0,'current_submission_eligible_count':0,
 'candidate_ipSAE_available_at_update':False}
ledger['bc2_v1_score_subset']='One rejected candidate for scoring-adapter diagnostic only; all10 manifest retained remotely'
ledger['old_dependency_jobs_cancelled_for_recovery']=['64536757','64536758','64536761','64536762','64536764','64536765']
ledger['dependency_policy']='Success-only with kill-on-invalid-dependency for new jobs; BC2 v1 diagnostic runs independently; at most two candidate-scoring workers'
path.write_text(json.dumps(ledger,indent=2))
print('RECOVERY_LEDGER_UPDATED')
