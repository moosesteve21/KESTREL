#!/usr/bin/env bash
# run_chain.sh - run several batches' searches back to back, pausing for review only when
# a clean signal turns up.
#
# Usage (in the tess tmux session):
#     ./run_chain.sh b037 b038 b039 ...
#
# For each batch in turn: ./pipeline.sh targets/<b>.txt results/<b> 8 "batch <b>" (the search;
# its follow-up then runs in its own fu-<b> window, as usual). When a search finishes:
#   - no clean signals: the next batch's search starts;
#   - one or more clean signals: the chain stops and posts to Discord (search channel, if configured)
#     that it's paused for review, with the command to resume.
# It also stops on any problem (pipeline.sh not exiting 0: a search incomplete or stopped,
# a failure, a typo guard), and before the next batch if results/STOP_CHAIN exists (the file
# is then removed). The search alone decides "clean"; later checks only downgrade signals.
# The user's rule (2026-10-15): batches may run as a chain only while no clean signal appears.
# Settings: CHAIN_WORKERS (default 8); CHAIN_TARGETS / CHAIN_RESULTS (folders, for tests).
set -uo pipefail

[[ $# -ge 1 ]] || { echo "usage: $0 batch [batch ...]" >&2; exit 1; }
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
tdir=${CHAIN_TARGETS:-targets}; rdir=${CHAIN_RESULTS:-results}
stop="$rdir/STOP_CHAIN"
notify() { "$here/server/notify.sh" --kind search "$*"; echo "run_chain: $*"; }
batches=("$@")

for i in "${!batches[@]}"; do
    b=${batches[$i]}
    rest=("${batches[@]:$i}")
    if [[ -f "$stop" ]]; then
        rm -f "$stop"
        notify "chain stopped by $stop before $b; resume with ./run_chain.sh ${rest[*]}"
        exit 0
    fi
    "$here/pipeline.sh" "$tdir/$b.txt" "$rdir/$b" "${CHAIN_WORKERS:-8}" "batch $b"
    rc=$?
    if (( rc != 0 )); then
        notify "chain stopped: $b ended with exit $rc (see $rdir/$b/pipeline.log); resume with ./run_chain.sh ${rest[*]}"
        exit "$rc"
    fi
    clean=$("$here/.venv/bin/python" -c "
import csv, sys
rows = csv.DictReader(open(sys.argv[1]))
print(sum(1 for r in rows if r.get('status') == 'clean'))" "$rdir/$b/candidates_all.csv" 2>/dev/null || echo 0)
    if (( clean > 0 )); then
        next=("${batches[@]:$((i + 1))}")
        notify "$b: $clean clean signal(s): chain paused for review$( (( ${#next[@]} )) && echo "; resume with ./run_chain.sh ${next[*]}")"
        exit 0
    fi
done
notify "chain finished: ${batches[*]} (no clean signals; follow-ups may still be running)"
