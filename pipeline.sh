#!/usr/bin/env bash
# pipeline.sh - search a target list and run every follow-up step, unattended.
#
# Usage:
#     ./pipeline.sh targets/<name>.txt results/<name> [n_workers] ["note for the database"]
#
# Steps, each logged to <out>/pipeline.log; the pipeline stops at the first failure:
#   1. run_parallel.sh       search (retries unfinished stars; must end complete)
#   2. enrich_candidates.py  catalogue checks           } skipped when the search
#   3. pixel_check.py        pixel-level checks         } logged no signals
#   4. reduction_check.py    independent-reduction check}
#   4a. vsx_shortlist.py     variable stars for VSX (uses the pixel and reduction verdicts)
#   4b. transit_fit.py       physical fit, density check } clean candidates; need the
#   4c. fpp_check.py         TRICERATOPS FPP             } vetting env (.venv-vet)
#   5. tessdb.py ingest      load the run into db/tess.db (also every 100 stars
#                            during the search: rolling updates, ROLLING_INGEST)
#   6. tessdb.py backup      dated copy in db/backups/, and records/*.csv (notes,
#                            watchlist, batches: the tables only the database holds)
# Rerunning the same command resumes: the search skips finished stars and the
# later steps are cheap to repeat (ingest replaces the run's rows).
#
# Exit status: 0 done; 3 search incomplete after its retries (rerun to retry);
# 4 the folder is busy; 6 the search was stopped (<out>/STOP or interrupted; rerun
# to continue); otherwise the exit status of the step that failed.
# To pause a long search: touch <out>/STOP. Workers finish their current star.
set -uo pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: $0 targets_file out_dir [n_workers] [note]" >&2
    exit 1
fi
targets=$1; root=$2; n=${3:-8}; note=${4:-"pipeline run on $(basename "$targets")"}
here=$(cd "$(dirname "$0")" && pwd)
py="$here/.venv/bin/python"
[[ -x "$py" ]] || py=python
cd "$here"

mkdir -p "$root"
exec 8> "$root/.pipeline.lock"
if ! flock -n 8; then
    echo "another pipeline is already running on $root" >&2
    exit 4
fi
logf="$root/pipeline.log"
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$logf"; }

step() {                     # step NAME COMMAND...
    local name=$1; shift
    log "== $name: $*"
    "$@" >> "$logf" 2>&1
    local rc=$?
    if (( rc != 0 )); then
        log "!! $name failed (exit $rc); see $logf"
        exit "$rc"
    fi
    log "   $name done"
}

log "pipeline start: $targets -> $root, $n workers"

# when the pipeline ends, however it ends, post the outcome to Discord (server/notify.sh;
# does nothing unless ~/.config/planetfinder/discord exists)
finish() {
    local rc=$? name msg db="${TESSDB:-$here/db/tess.db}"
    name=$(basename "$root")
    case $rc in
        0) msg="$name complete: $(sqlite3 "$db" "SELECT r.n_targets || ' stars, ' ||
               (SELECT count(*) FROM candidates c WHERE c.run_id = r.run_id AND c.status != 'same_period') || ' signals, ' ||
               (SELECT count(*) FROM candidates c WHERE c.run_id = r.run_id AND c.status = 'clean') || ' clean, ' ||
               (SELECT count(*) FROM vsx_shortlist v WHERE v.run_id = r.run_id) || ' for VSX'
               FROM runs r WHERE r.name = '$name'" 2>/dev/null)" ;;
        6) msg="$name stopped: $(wc -l < "$root/failed_targets.txt" 2>/dev/null || echo ?) star(s) left; rerun to continue" ;;
        3) msg="$name search incomplete: $(wc -l < "$root/failed_targets.txt" 2>/dev/null || echo ?) star(s) failed; rerun to retry" ;;
        *) msg="$name FAILED: $(grep ' !! ' "$logf" | tail -1 | sed 's/.* !! //' | cut -c1-200)" ;;
    esac
    "$here/server/notify.sh" "$msg (after $((SECONDS / 3600))h $((SECONDS % 3600 / 60))m)"
}
trap finish EXIT

# 1. search: run_parallel.sh retries by itself; anything but 0 stops the pipeline.
# Rolling updates: every ROLLING_INGEST stars (default 100; 0 turns them off) the
# results so far are loaded into the database, so a long run can be browsed as it goes.
export ROLLING_INGEST=${ROLLING_INGEST:-100} RUN_NOTE="$note"
log "== search: run_parallel.sh $targets $root $n"
"$here/run_parallel.sh" "$targets" "$root" "$n" 2>&1 | tee -a "$logf"
rc=${PIPESTATUS[0]}
if (( rc == 3 )); then
    log "!! search incomplete after retries: $(wc -l < "$root/failed_targets.txt") star(s) in $root/failed_targets.txt"
    log "   rerun this command to retry them; the follow-up steps wait until the search is complete"
    exit 3
elif (( rc == 6 )); then
    log "!! search stopped (STOP file or interrupted): $(wc -l < "$root/failed_targets.txt") star(s) left"
    log "   rerun this command to continue; the follow-up steps wait until the search is complete"
    exit 6
elif (( rc != 0 )); then
    log "!! search failed (exit $rc)"
    exit "$rc"
fi
[[ -f "$root/attention.txt" ]] && log "   note: $(wc -l < "$root/attention.txt") star(s) need attention (search stopped partway): $root/attention.txt"

# 2-4. follow-up checks, only if there is something to check
signals=$(( $(wc -l < "$root/candidates_all.csv") - 1 ))
if (( signals > 0 )); then
    step enrich "$py" enrich_candidates.py "$root/candidates_all.csv"
    step pixel_check "$py" pixel_check.py "$root/enriched.csv"
    step reduction_check "$py" reduction_check.py "$root/enriched.csv"
    step vsx_shortlist "$py" vsx_shortlist.py --no-db "$root"
    if [[ -x "$here/.venv-vet/bin/python" ]]; then
        step transit_fit "$py" transit_fit.py "$root/enriched.csv"
        step fpp_check "$py" fpp_check.py "$root/enriched.csv"
    else
        log "   no vetting environment (.venv-vet): transit fit and FPP check skipped"
    fi
else
    log "   no signals logged: catalogue, pixel and reduction checks skipped"
fi

# 5-6. database
step ingest "$py" tessdb.py ingest "$root" --note "$note"
step backup "$py" tessdb.py backup
log "pipeline complete: $root"
"$py" tessdb.py list 2>/dev/null | grep -E "^run|^$(basename "$root") " | tee -a "$logf"
exit 0
