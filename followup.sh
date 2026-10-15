#!/usr/bin/env bash
# followup.sh - every step after the search, for one batch: catalogue, pixel and reduction
# checks, VSX shortlist, transit fit, depth map, FPP, then the database ingest and backup.
#
# Usage:
#     ./followup.sh results/<name> ["note for the database"]
#
# pipeline.sh starts this by itself when a search finishes (in a tess tmux window named
# fu-<name>), so the next batch's search can start while this one, mostly waiting on MAST,
# carries on. Run it by hand to redo a batch's follow-up (the steps are cheap to repeat;
# ingest replaces the run's rows).
#   - One follow-up at a time: a later batch waits for the lock (followup.lock beside the
#     database) and goes when the earlier one is done.
#   - The batch folder's lock (<out>/.pipeline.lock) is held too, so a rerun of pipeline.sh
#     on the same folder can't overlap it.
#   - Its steps run at low CPU priority (FOLLOWUP_NICE, default 10), so a search running
#     alongside keeps its speed; the transit fit and TRICERATOPS are the CPU-heavy ones.
# Posts the outcome to Discord (server/notify.sh), as pipeline.sh used to.
set -uo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 out_dir [note]" >&2
    exit 1
fi
root=$1; note=${2:-"pipeline run on $(basename "$root")"}
here=$(cd "$(dirname "$0")" && pwd)
py="$here/.venv/bin/python"
[[ -x "$py" ]] || py=python
cd "$here"
db="${TESSDB:-$here/db/tess.db}"
logf="$root/pipeline.log"
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$logf"; }

exec 7> "$(dirname "$db")/followup.lock"
if ! flock -n 7; then
    log "follow-up of $root waiting for another batch's follow-up to finish"
    flock 7
fi
exec 8> "$root/.pipeline.lock"
if ! flock -n 8; then
    log "follow-up of $root waiting for the folder (a pipeline.sh still running on it)"
    flock 8
fi
SECONDS=0
log "follow-up start: $root"

step() {                     # step NAME COMMAND...
    local name=$1; shift
    log "== $name: $*"
    nice -n "${FOLLOWUP_NICE:-10}" "$@" >> "$logf" 2>&1
    local rc=$?
    if (( rc != 0 )); then
        log "!! $name failed (exit $rc); see $logf"
        exit "$rc"
    fi
    log "   $name done"
}

finish() {
    local rc=$? name msg
    name=$(basename "$root")
    if (( rc == 0 )); then
        msg="$name complete: $(sqlite3 "$db" "SELECT r.n_targets || ' stars, ' ||
               (SELECT count(*) FROM candidates c WHERE c.run_id = r.run_id AND c.status != 'same_period') || ' signals, ' ||
               (SELECT count(*) FROM candidates c WHERE c.run_id = r.run_id AND c.status = 'clean') || ' clean, ' ||
               (SELECT count(*) FROM vsx_shortlist v WHERE v.run_id = r.run_id) || ' for VSX'
               FROM runs r WHERE r.name = '$name'" 2>/dev/null)"
    else
        msg="$name follow-up FAILED: $(grep ' !! ' "$logf" | tail -1 | sed 's/.* !! //' | cut -c1-200)"
    fi
    "$here/server/notify.sh" --kind batch "$msg (follow-up $((SECONDS / 3600))h $((SECONDS % 3600 / 60))m)"
}
trap finish EXIT

signals=$(( $(wc -l < "$root/candidates_all.csv") - 1 ))
if (( signals > 0 )); then
    step enrich "$py" enrich_candidates.py "$root/candidates_all.csv"
    step pixel_check "$py" pixel_check.py "$root/enriched.csv"
    step reduction_check "$py" reduction_check.py "$root/enriched.csv"
    step vsx_shortlist "$py" vsx_shortlist.py --no-db "$root"
    if [[ -x "$here/.venv-vet/bin/python" ]]; then
        step transit_fit "$py" transit_fit.py "$root/enriched.csv"
        step depth_map "$py" depth_map_check.py "$root/enriched.csv"
        step fpp_check "$py" fpp_check.py "$root/enriched.csv"
    else
        log "   no vetting environment (.venv-vet): transit fit, depth map and FPP check skipped"
    fi
else
    log "   no signals logged: catalogue, pixel and reduction checks skipped"
fi

step ingest "$py" tessdb.py ingest "$root" --note "$note"
step backup "$py" tessdb.py backup
log "pipeline complete: $root"
"$py" tessdb.py list 2>/dev/null | grep -E "^run|^$(basename "$root") " | tee -a "$logf"
exit 0
