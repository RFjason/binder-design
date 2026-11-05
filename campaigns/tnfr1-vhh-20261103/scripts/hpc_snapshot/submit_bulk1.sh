#!/bin/bash
# Login-node management only: Slurm submission, dependency setup and file bookkeeping.
set -euo pipefail
ROOT=/dssg/home/acct-clshzh/jasoncao1280/benchmarks/tnfr1_vhh_20261103
grep -q '"ok": true' "$ROOT/reports/bulk1_preflight.json"
IDS="$ROOT/reports/bulk1_submission_ids.tsv"
if [[ "${1:-}" == '--resume-held' ]]; then
    while IFS=$'\t' read -r stage job; do
        case "$stage" in
            generation) GEN=$job ;;
            collection) COLLECT=$job ;;
            scoring) SCORE=$job ;;
            audit) AUDIT=$job ;;
        esac
    done < "$IDS"
    for job in "$GEN" "$COLLECT" "$SCORE" "$AUDIT"; do
        scontrol show job "$job" | grep -q 'Reason=JobHeldUser'
    done
elif ! mkdir "$ROOT/reports/bulk1_submission.lock"; then
    echo 'Bulk1 submission lock already exists; inspect the existing graph before retrying.' >&2
    exit 2
fi
trap 'echo "Submission setup interrupted. Inspect held jobs and reports/bulk1_submission_ids.tsv before recovery." >&2' ERR
if [[ "${1:-}" != '--resume-held' ]]; then
: > "$IDS"
GEN=$(sbatch --parsable --hold "$ROOT/scripts/bulk1_generation.slurm")
printf 'generation\t%s\n' "$GEN" >> "$IDS"
COLLECT=$(sbatch --parsable --hold "$ROOT/scripts/bulk1_collection.slurm")
printf 'collection\t%s\n' "$COLLECT" >> "$IDS"
SCORE=$(sbatch --parsable --hold --kill-on-invalid-dep=yes "$ROOT/scripts/bulk1_scoring.slurm")
printf 'scoring\t%s\n' "$SCORE" >> "$IDS"
AUDIT=$(sbatch --parsable --hold "$ROOT/scripts/bulk1_audit.slurm")
printf 'audit\t%s\n' "$AUDIT" >> "$IDS"
fi
while IFS=$'\t' read -r task method ident genpart scorepart; do
    scontrol update "JobId=${GEN}_$task" "Partition=$genpart"
    scontrol update "JobId=${COLLECT}_$task" "Dependency=afterany:${GEN}_$task"
    scontrol update "JobId=${SCORE}_$task" "Partition=$scorepart" "Dependency=afterok:${COLLECT}_$task"
    scontrol update "JobId=${AUDIT}_$task" "Dependency=afterany:${SCORE}_$task"
done < "$ROOT/inputs/bulk1/batches.tsv"
SUMMARY=$(sbatch --parsable --dependency="afterany:$AUDIT" "$ROOT/scripts/bulk1_summary.slurm")
printf 'summary\t%s\n' "$SUMMARY" >> "$IDS"
printf '{"campaign":"bulk1","generation":"%s","collection":"%s","scoring":"%s","audit":"%s","summary":"%s","generation_throttle":6,"scoring_throttle":2,"generation_walltime":"02:00:00","scoring_walltime":"01:00:00"}\n' "$GEN" "$COLLECT" "$SCORE" "$AUDIT" "$SUMMARY" > "$ROOT/reports/bulk1_submission.json"
scontrol release "$AUDIT"
scontrol release "$SCORE"
scontrol release "$COLLECT"
scontrol release "$GEN"
trap - ERR
cat "$ROOT/reports/bulk1_submission.json"
