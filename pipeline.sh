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
#   4c. depth_map_check.py   per-source pixel depth map  } vetting env (.venv-vet);
#   4d. fpp_check.py         TRICERATOPS FPP             } the depth map uses the fit's
#                            ephemeris, and the FPP skips dips it puts off target
#   5. tessdb.py ingest      load the run into db/tess.db (also every 100 stars
#                            during the search: rolling updates, ROLLING_INGEST)
#   6. tessdb.py backup      dated copy in db/backups/, and records/*.csv (notes,
#                            watchlist, batches: the tables only the database holds)
# Steps 2-6 are followup.sh. When the search is complete this script posts "search done"
# and hands the batch to followup.sh in a tess tmux window (fu-<name>), so the next
# batch's search can start while this one's network-bound follow-up runs; follow-ups
# queue and go one at a time, at low CPU priority. FOLLOWUP=inline runs them here
# instead (as before; used by the tests).
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

# Guard against typos (e.g. targets/b029.txt with results/b028), and against running over a
# batch that's finished or whose follow-up is still pending. FORCE=1 skips the last two;
# ANY_NAMES=1 allows names that differ (scratch tests).
fail() { echo "pipeline.sh: $*" >&2; exit 2; }
tname=$(basename "$targets"); tname=${tname%.*}
[[ -f "$targets" ]] || fail "no targets file $targets"
if [[ "$tname" != "$(basename "$root")" && -z "${ANY_NAMES:-}" ]]; then
    fail "targets $targets don't match results folder $root (expected results/$tname); ANY_NAMES=1 to allow"
fi
if [[ -f "$root/all_targets.txt" ]]; then
    want=$(sed 's/#.*//' "$targets" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | grep -v '^$' | sort -u)
    have=$(sort -u "$root/all_targets.txt")
    [[ "$want" == "$have" ]] || fail "$root already holds a run of different stars ($(wc -l < "$root/all_targets.txt") targets): not overwritten"
fi
if [[ -f "$root/pipeline.log" && -z "${FORCE:-}" ]]; then
    last=$(grep -E 'pipeline complete:|search complete: follow-up queued' "$root/pipeline.log" | tail -1)
    case $last in
        *"pipeline complete:"*) fail "$root is already complete ($last); FORCE=1 to run it again" ;;
        *"search complete: follow-up queued"*) fail "$root: search done and its follow-up is queued or running (tmux window fu-$(basename "$root")); FORCE=1 to run it again" ;;
    esac
fi

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
        0) [[ -n "${handed_over:-}" ]] || return 0      # inline follow-up posts its own message
           msg="$name search done: $(grep -vc ',same_period,' "$root/candidates_all.csv" 2>/dev/null | awk '{print $1 - 1}') signals; follow-up $handed_over; next batch can start" ;;
        6) msg="$name stopped: $(wc -l < "$root/failed_targets.txt" 2>/dev/null || echo ?) star(s) left; rerun to continue" ;;
        3) msg="$name search incomplete: $(wc -l < "$root/failed_targets.txt" 2>/dev/null || echo ?) star(s) failed; rerun to retry" ;;
        *) msg="$name FAILED: $(grep ' !! ' "$logf" | tail -1 | sed 's/.* !! //' | cut -c1-200)" ;;
    esac
    "$here/server/notify.sh" --kind search "$msg (after $((SECONDS / 3600))h $((SECONDS % 3600 / 60))m)"
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

# 2-6. follow-up: in its own tess window (default), or here (FOLLOWUP=inline)
exec 8>&-                                     # release the folder lock: followup.sh takes it
if [[ "${FOLLOWUP:-queue}" == inline ]]; then
    trap - EXIT
    "$here/followup.sh" "$root" "$note"
    exit $?
fi
fu=(env)                                      # carry the settings that matter (tmux starts a fresh shell)
for v in TESSDB PLANETFINDER_DISCORD NOTIFY_DRYRUN FOLLOWUP_NICE PLANETFINDER_LOG; do
    [[ -n "${!v:-}" ]] && fu+=("$v=${!v}")
done
fu+=("$here/followup.sh" "$root" "$note")
if tmux has-session -t "=tess" 2>/dev/null && [[ "${FOLLOWUP_TMUX:-1}" == 1 ]]; then
    tmux new-window -d -t "=tess:" -n "fu-$(basename "$root")" -c "$here" "$(printf '%q ' "${fu[@]}")"
    handed_over="queued (tess window fu-$(basename "$root"))"
else
    setsid nohup "${fu[@]}" > /dev/null 2>&1 < /dev/null &
    handed_over="queued (background process $!)"
fi
log "search complete: follow-up $handed_over; its steps are logged here"
exit 0
