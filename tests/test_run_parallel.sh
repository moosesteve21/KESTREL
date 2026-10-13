#!/usr/bin/env bash
# tests/test_run_parallel.sh - check run_parallel's queue, retry and completeness
# logic in seconds, using tests/fake_search.py in place of find_planets.py.
# Run from anywhere.
set -u
here=$(cd "$(dirname "$0")/.." && pwd)
cd "$here"
py="$here/.venv/bin/python"; [[ -x "$py" ]] || py=python3
tmp=$(mktemp -d)
cp tests/fake_search.py ./_fake_search_test.py
trap 'rm -rf "$tmp" ./_fake_search_test.py' EXIT
export STAGGER=0 POLL=0.05 RETRY_WAIT=0
fails=0
pass() { echo "ok   $1"; }
fail() { echo "FAIL $1"; fails=$((fails + 1)); }
run() {        # run NAME WORKERS -> exit status of run_parallel.sh on $tmp/NAME.txt
    ./run_parallel.sh "$tmp/$1.txt" "$tmp/$1" "$2" _fake_search_test.py >> "$tmp/$1.log" 2>&1
}
check() {      # check NAME EXPECTED_EXIT WORKERS TARGETS...
    local name=$1 want=$2 n=$3; shift 3
    printf "%s\n" "$@" > "$tmp/$name.txt"
    run "$name" "$n"; local got=$?
    [[ $got == "$want" ]] && pass "$name (exit $got)" || fail "$name: exit $got, expected $want (log: $tmp/$name.log)"
}
attempts() {   # attempts DIR STAR -> attempts recorded across all workers
    cat "$1"/w*/attempts.txt 2>/dev/null | awk -F'\t' -v s="$2" '$1 == s' | wc -l
}
rows() {       # rows DIR STAR -> rows for STAR in candidates_all.csv
    awk -F, -v s="$2" '$1 == s' "$1/candidates_all.csv" | wc -l
}

# --- retries, failures, partial results, crashes
RETRY_WAIT=1 check flaky_recovers 0 2 STAR_1 FLAKY_2 STAR_3
[[ $(attempts "$tmp/flaky_recovers" FLAKY_2) == 2 ]] && pass "flaky star retried once" || fail "flaky star attempts"
grep -q "^FLAKY_2,1,attempt2" "$tmp/flaky_recovers/candidates_all.csv" && [[ $(rows "$tmp/flaky_recovers" FLAKY_2) == 1 ]] \
    && pass "only the successful attempt's rows are merged" || fail "flaky star rows in candidates_all.csv"
MAX_ATTEMPTS=2 check broken_fails 3 2 STAR_1 BROKEN_2
grep -qx BROKEN_2 "$tmp/broken_fails/failed_targets.txt" && [[ $(attempts "$tmp/broken_fails" BROKEN_2) == 2 ]] \
    && pass "broken star failed after 2 attempts" || fail "broken star not failed correctly"
check partial_attention 0 2 STAR_1 PARTIAL_2
grep -q PARTIAL_2 "$tmp/partial_attention/attention.txt" && pass "partial star in attention.txt" || fail "partial star not in attention.txt"
check crash_retried 0 2 CRASH_1 STAR_2 STAR_3
[[ $(attempts "$tmp/crash_retried" CRASH_1) == 2 && $(rows "$tmp/crash_retried" CRASH_1) == 1 ]] \
    && grep -q "^CRASH_1,1,attempt2" "$tmp/crash_retried/candidates_all.csv" \
    && pass "crashed worker's star re-queued and merged once" || fail "crashed star handling"

# --- a second run on a busy folder is refused
mkdir -p "$tmp/busy"; exec 7> "$tmp/busy/.run.lock"; flock -n 7
printf "STAR_1\n" > "$tmp/busy.txt"
run busy 2; [[ $? == 4 ]] && pass "busy_folder (exit 4)" || fail "busy folder not refused"
exec 7>&-

# --- many workers, many stars: every star searched exactly once
seq -f "STAR_%04g" 1 500 > "$tmp/many.txt"
FAKE_SLEEP=0.01 run many 8; rc=$?
dupes=$(cat "$tmp"/many/w*/attempts.txt | cut -f1 | sort | uniq -d | wc -l)
total=$(cat "$tmp"/many/w*/attempts.txt | wc -l)
used=$(ls -d "$tmp"/many/w*/ | wc -l)
[[ $rc == 0 && $dupes == 0 && $total == 500 && $(wc -l < "$tmp/many/done_all.txt") == 500 \
   && $(($(wc -l < "$tmp/many/candidates_all.csv") - 1)) == 500 && $used == 8 ]] \
    && pass "500 stars, 8 workers: each searched once, all merged" \
    || fail "concurrency: exit $rc, $total attempts, $dupes duplicated, $used workers"

# --- STOP pauses the run; rerunning continues without repeating stars
printf "STAR_1\nSTAR_2\nSTOPPER_3\nSTAR_4\nSTAR_5\nSTAR_6\n" > "$tmp/stop.txt"
run stop 1; rc=$?
left=$(wc -l < "$tmp/stop/failed_targets.txt" 2>/dev/null || echo 0)
[[ $rc == 6 && $left == 3 ]] && pass "STOP file pauses the run (exit 6, 3 left)" || fail "stop: exit $rc, $left left"
run stop 2; rc=$?
[[ $rc == 0 && ! -f "$tmp/stop/STOP" && $(cat "$tmp"/stop/w*/attempts.txt | cut -f1 | sort | uniq -d | wc -l) == 0 ]] \
    && pass "rerun after STOP completes, nothing searched twice, worker count changed" \
    || fail "resume after stop: exit $rc"

# --- a claim left by a dead worker (supervisor killed) is returned at the next start
printf "STAR_1\nSTAR_2\n" > "$tmp/dead.txt"
mkdir -p "$tmp/dead"
"$py" -c "
import sys; sys.path.insert(0, '.')
from workqueue import Queue
q = Queue('$tmp/dead/queue.db'); q.add(['STAR_1', 'STAR_2'])
q.claim(0, pid=999999)          # no such process
"
run dead 2; rc=$?
[[ $rc == 0 && $(attempts "$tmp/dead" STAR_1) == 1 ]] && grep -q "claim of a dead worker" "$tmp/dead.log" \
    && pass "dead worker's claim returned and searched" || fail "dead claim: exit $rc"

# --- a claim whose pid now belongs to an unrelated live process (pid reuse after a
#     reboot) is returned too, instead of being waited on forever
printf "STAR_1\nSTAR_2\n" > "$tmp/reused.txt"
mkdir -p "$tmp/reused"
sleep 60 & stranger=$!
"$py" -c "
import sys; sys.path.insert(0, '.')
from workqueue import Queue
q = Queue('$tmp/reused/queue.db'); q.add(['STAR_1', 'STAR_2'])
q.claim(0, pid=$stranger)       # alive, but not one of our workers
"
run reused 2; rc=$?
kill $stranger 2>/dev/null
[[ $rc == 0 && $(attempts "$tmp/reused" STAR_1) == 1 ]] && grep -q "claim of a dead worker" "$tmp/reused.log" \
    && pass "claim held by an unrelated live process returned and searched" || fail "reused pid claim: exit $rc"

# --- a folder from the old round-robin runner is resumed from its done.txt files
mkdir -p "$tmp/old/w0"
printf "STAR_1\tok_0_candidates\t2026-01-01T00:00:00\nSTAR_2\terror: X\t2026-01-01T00:00:00\n" > "$tmp/old/w0/done.txt"
printf "STAR_1\nSTAR_2\nSTAR_3\n" > "$tmp/old.txt"
run old 2; rc=$?
[[ $rc == 0 && $(attempts "$tmp/old" STAR_1) == 0 && $(attempts "$tmp/old" STAR_2) == 1 ]] \
    && pass "old-style folder: finished stars kept, unfinished searched" || fail "old folder: exit $rc"

# --- a rerun gives failed stars fresh attempts
MAX_ATTEMPTS=2 run broken_fails 2; rc=$?
[[ $rc == 3 && $(attempts "$tmp/broken_fails" BROKEN_2) == 4 ]] && pass "rerun retries failed stars" \
    || fail "rerun of failed stars: exit $rc, $(attempts "$tmp/broken_fails" BROKEN_2) attempts"

# --- rolling updates: the results so far go into the database every N stars
seq -f "ROLL_%03g" 1 120 > "$tmp/roll.txt"
TESSDB="$tmp/roll.db" ROLLING_INGEST=50 RUN_NOTE="rolling test" FAKE_SLEEP=0.02 run roll 4; rc=$?
n_roll=$(grep -c "rolling update: " "$tmp/roll.log")
n_db=$("$py" -c "
import sqlite3; c = sqlite3.connect('$tmp/roll.db')
r = c.execute('SELECT note, n_targets FROM runs').fetchone()
print(f'{r[1]}|{r[0]}')" 2>/dev/null)
n_at=${n_db%%|*}       # stars loaded at the last update: 100-119 (several can finish between checks)
[[ $rc == 0 && $n_roll == 2 && $n_at -ge 100 && $n_at -lt 120 && $n_db == *"rolling test [in progress: $n_at/120 searched]" ]] \
    && pass "rolling updates every 50 stars (2 loads; database at $n_at/120)" \
    || fail "rolling updates: exit $rc, $n_roll loads, database '$n_db'"

# --- a stopped run is loaded into the database on the way out
printf "STAR_1\nSTAR_2\nSTOPPER_3\nSTAR_4\nSTAR_5\n" > "$tmp/stopdb.txt"
TESSDB="$tmp/stopdb.db" ROLLING_INGEST=100 RUN_NOTE="stop test" run stopdb 1; rc=$?
note=$("$py" -c "
import sqlite3; print(sqlite3.connect('$tmp/stopdb.db').execute('SELECT note FROM runs').fetchone()[0])" 2>/dev/null)
[[ $rc == 6 && $note == "stop test [stopped: 3/5 searched]" ]] && pass "stopped run loaded into the database ($note)" \
    || fail "stopped run: exit $rc, note '$note'"

# --- a search started while prefetch.py holds the folder's lock waits for it
mkdir -p "$tmp/pfwait"; printf "STAR_1\nSTAR_2\n" > "$tmp/pfwait.txt"
( exec 9> "$tmp/pfwait/.prefetch.lock"; flock -x 9; sleep 2 ) &
sleep 0.3
PREFETCH_POLL=0.2 run pfwait 2; rc=$?
[[ $rc == 0 ]] && grep -q "waiting for the prefetch" "$tmp/pfwait.log" && grep -q "prefetch finished" "$tmp/pfwait.log" \
    && pass "search waits for a running prefetch, then runs" || fail "prefetch wait: exit $rc"

# --- a second prefetch on the same folder is refused (before any MAST query)
mkdir -p "$tmp/pfbusy"; printf "STAR_1\n" > "$tmp/pfbusy.txt"
( exec 9> "$tmp/pfbusy/.prefetch.lock"; flock -x 9; sleep 20 ) & holder=$!
sleep 0.3
"$py" prefetch.py "$tmp/pfbusy.txt" "$tmp/pfbusy" >> "$tmp/pfbusy.log" 2>&1; rc=$?
kill $holder 2>/dev/null
[[ $rc == 4 ]] && grep -q "already running" "$tmp/pfbusy.log" && pass "second prefetch on a folder refused (exit 4)" \
    || fail "second prefetch: exit $rc"

# --- a manifest line still being written is skipped, not an error
printf '%s\n%s' '{"star": "STAR_A", "status": "no_data", "files": []}' '{"star": "STAR_B", "sta' > "$tmp/half.jsonl"
out=$(PREFETCH_MANIFEST="$tmp/half.jsonl" "$py" -c "
import sys; sys.path.insert(0, '.'); import find_planets as fp
print(fp.prefetched('STAR_A')['status'], fp.prefetched('STAR_B'))" 2>/dev/null | tail -1)
[[ $out == "no_data None" ]] && pass "half-written manifest line skipped" || fail "half-written manifest line: '$out'"

(( fails == 0 )) && echo "all tests passed" || echo "$fails test(s) failed"
exit $(( fails > 0 ))
