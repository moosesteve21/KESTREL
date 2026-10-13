"""Stand-in for find_planets.py to test run_parallel's queue, retry and completeness
logic quickly and deterministically. Worker mode only (--queue Q --worker K --out D),
like find_planets.py under run_parallel.py. Behaviour by target name:
  FLAKY_*   error on the first attempt, ok afterwards
  BROKEN_*  always error
  PARTIAL_* partial_..._then_error (counts as finished, needs attention)
  CRASH_*   the first attempt kills the worker (exit 1) after claiming the star,
            leaving it marked running
  STOPPER_* creates the run's STOP file, then finishes ok
  anything else: ok
Every star gets one row in candidates.csv per attempt (so leftover rows from a
failed attempt can be checked), and every attempt is appended to <out>/attempts.txt
as "star<TAB>worker<TAB>attempt". FAKE_SLEEP (seconds) simulates search time."""
import argparse
import csv
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from workqueue import Queue  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--queue", required=True)
ap.add_argument("--worker", type=int, required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
out = Path(args.out)
out.mkdir(parents=True, exist_ok=True)
q = Queue(args.queue)
stop = Path(args.queue).parent / "STOP"
sleep = float(os.environ.get("FAKE_SLEEP", 0))
while not stop.exists():
    t = q.claim(args.worker, os.getpid())
    if t is None:
        break
    n = q.con.execute("SELECT attempts FROM queue WHERE star = ?", (t,)).fetchone()[0]
    cf = out / "candidates.csv"                 # like find_planets: drop this worker's
    if cf.exists():                             # rows from an earlier attempt
        keep = [r for r in csv.reader(open(cf)) if r and r[0] != t]
        with open(cf, "w", newline="") as f:
            csv.writer(f).writerows(keep)
    with open(out / "attempts.txt", "a") as f:
        f.write(f"{t}\t{args.worker}\t{n}\n")
    time.sleep(sleep)
    new = not (out / "candidates.csv").exists() or not (out / "candidates.csv").stat().st_size
    with open(out / "candidates.csv", "a", newline="") as f:
        w = csv.writer(f)
        if new:                                 # the columns tessdb.py ingest reads
            w.writerow(["target", "candidate", "status", "period_d", "t0_btjd", "duration_h",
                        "depth_ppm", "snr", "sde", "flags", "plot"])
        w.writerow([t, 1, f"attempt{n}", 1.0, 0.0, 2.0, 100, 10, 10, "", ""])
    if t.startswith("CRASH_") and n == 1:
        os._exit(1)                         # dies holding the claim
    if t.startswith("STOPPER_"):
        stop.touch()
    if t.startswith("BROKEN_") or (t.startswith("FLAKY_") and n == 1):
        status = "error: SimulatedError: test"
    elif t.startswith("PARTIAL_"):
        status = "partial_1_candidates_then_error"
    else:
        status = "ok_1_candidates"
    with open(out / "done.txt", "a") as f:
        f.write(f"{t}\t{status}\t{time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
    q.finish(t, status)
sys.exit(0)
