#!/usr/bin/env bash
# server/batch_status.sh - list batches whose pipeline didn't finish (stopped, failed,
# or cut off by a reboot), with the command that resumes each one.
cd "$(dirname "$0")/.." || exit 1
found=0
for log in results/*/pipeline.log; do
    [[ -f "$log" ]] || continue
    root=$(dirname "$log"); name=$(basename "$root")
    start=$(grep "pipeline start:" "$log" | tail -1)
    [[ -z "$start" ]] && continue
    # finished if the last start is followed by "pipeline complete"
    tail -n +"$(grep -n "pipeline start:" "$log" | tail -1 | cut -d: -f1)" "$log" \
        | grep -q "pipeline complete" && continue
    # still running?
    pgrep -f "run_parallel.py [^ ]* $root " >/dev/null && { echo "$name: running now"; continue; }
    targets=$(sed -E 's/.*pipeline start: ([^ ]+) -> .*/\1/' <<< "$start")
    workers=$(sed -E 's/.*, ([0-9]+) workers.*/\1/' <<< "$start")
    done_n=""
    [[ -f "$root/queue.db" ]] && done_n=$(sqlite3 "$root/queue.db" \
        "SELECT sum(state='done') || '/' || count(*) FROM queue" 2>/dev/null)
    echo "$name: did not finish (${done_n:-?} searched; last: $(tail -1 "$log" | cut -c1-100))"
    echo "    resume: ./pipeline.sh $targets $root $workers \"batch $name\""
    found=1
done
(( found )) || echo "nothing to resume"
