#!/usr/bin/env python3
"""
run_parallel.py - search a target list with several workers sharing one queue,
retrying stars that end in an error, until the run is complete.

Usage:
    ./run_parallel.sh targets/<name>.txt results/<name> [n_workers] [script]
    (run_parallel.sh calls this; defaults: 8 workers, find_planets.py)
    Environment: MAX_ATTEMPTS (default 3; MAX_PASSES is accepted too), RETRY_WAIT
    (seconds before an errored star is retried, default 300), STAGGER (seconds
    between worker starts, default 20).

The targets go into <out>/queue.db in file order (score order for batches).
Each worker claims the best-ranked star that is due, searches it, records the
outcome and claims the next, so no worker idles while stars are left. A star
that ends in an error is retried after RETRY_WAIT seconds, up to MAX_ATTEMPTS
attempts in all. A worker that dies has its star returned to the queue (that
counts as an attempt) and is replaced. Worker k writes only to <out>/w<k>.

Rerunning the same command resumes: finished stars are skipped, stars still
marked running by a dead worker are returned to the queue, and new stars in the
target file are added at the end. The worker count may differ between runs.
A folder from the older round-robin runner (w*/done.txt, no queue.db) is
resumed by seeding the queue from its done.txt files.

Pausing: create <out>/STOP. Workers finish their current star and exit, and the
run stops with status 6; rerunning the command removes STOP and continues.

Rolling updates: with ROLLING_INGEST=N (pipeline.sh sets 100), every N stars
finished the results so far are loaded into tess.db (tessdb.py ingest; the run's
note says "[in progress: x/y searched]"), so long runs can be browsed as they go;
a run stopped by STOP or a signal is loaded too ("[stopped: x/y searched]").
The final load, after the follow-up checks, replaces them.

At the end <out>/candidates_all.csv (each star's rows from the attempt that
finished it), done_all.txt (one line per star: star, status, time, worker),
failed_targets.txt and attention.txt are written from the queue.

Exit status:
  0  complete: every target has a final status
  3  incomplete: stars used all their attempts (listed in failed_targets.txt);
     running the same command again gives them fresh attempts
  4  another run is already using this output folder
  6  stopped (STOP file or interrupted) with stars left; rerun to continue
  1  bad arguments
<out>/attention.txt lists stars whose search stopped partway (partial_..._then_error):
they count as complete, but their candidate lists may be short.
"""
import csv
import datetime
import fcntl
import importlib.util
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from workqueue import Queue, is_error, now_iso  # noqa: E402

MAX_ATTEMPTS = int(os.environ.get("MAX_ATTEMPTS", os.environ.get("MAX_PASSES", 3)))
RETRY_WAIT = float(os.environ.get("RETRY_WAIT", 300))
STAGGER = float(os.environ.get("STAGGER", 20))
POLL = float(os.environ.get("POLL", 2))       # seconds between supervisor checks
PREFETCH_POLL = float(os.environ.get("PREFETCH_POLL", 10))   # checks while a prefetch runs
ROLLING_INGEST = int(os.environ.get("ROLLING_INGEST", 0))   # load into tess.db every N stars (0: off)
RUN_NOTE = os.environ.get("RUN_NOTE", "")       # the run's database note (pipeline.sh sets both)


def log(msg):
    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {msg}", flush=True)


def read_targets(path):
    out = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def git(*args):
    return subprocess.run(["git", "-C", str(HERE), *args], capture_output=True, text=True).stdout.strip()


def write_run_info(root, targets_file, script, n):
    """Record what produced this run (read by tessdb.py). The first launch writes
    it; a resume keeps it and appends its own commit, time and worker count."""
    path = root / "run_info.json"
    version = {"commit": git("rev-parse", "HEAD"),     # untracked files (e.g. new target lists) aren't changes
               "uncommitted_changes": bool(git("status", "--porcelain", "--untracked-files=no"))}
    if path.exists():
        info = json.loads(path.read_text())
        info.setdefault("resumes", []).append({"started": now_iso(), **version, "workers": n})
    else:
        try:
            spec = importlib.util.spec_from_file_location("search", HERE / script)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            settings = {k: getattr(mod, k) for k in dir(mod) if k.isupper() and not k.startswith("_")}
            settings = json.loads(json.dumps(settings, default=lambda o: getattr(o, "tolist", lambda: str(o))()))
        except (Exception, SystemExit) as e:          # a script may parse arguments on import
            settings = {"error": f"settings not captured: {type(e).__name__}"}
        info = {"started": now_iso(), **version, "targets_file": str(targets_file), "script": script,
                "workers": n, "runner": "queue", "max_attempts": MAX_ATTEMPTS, "retry_wait": RETRY_WAIT,
                "settings": settings}
    path.write_text(json.dumps(info, indent=1))
    return info["commit"] if "resumes" not in info else info["resumes"][-1]["commit"]


def seed_from_old_run(q, root):
    """A folder made by the round-robin runner: mark its finished stars done."""
    n = 0
    for f in sorted(root.glob("w*/done.txt")):
        worker = int(f.parent.name[1:]) if f.parent.name[1:].isdigit() else None
        for line in f.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2 and not is_error(parts[1]):
                q.mark_done(parts[0], parts[1], parts[2] if len(parts) > 2 else None, worker)
                n += 1
    return n


def pid_alive(pid):
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def our_worker(pid, queue_db):
    """True if pid is alive AND is a worker of this queue. After a reboot or a hard stop,
    a claim's pid may have been reused by an unrelated process: such a claim is dead
    too, or the batch would wait on it forever. Where /proc can't say (no access, not
    Linux), fall back to "alive"."""
    if not pid_alive(pid):
        return False
    try:
        args = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
        cwd = os.readlink(f"/proc/{pid}/cwd")
    except (OSError, ValueError):
        return True
    args = [a.decode(errors="replace") for a in args]
    if "--queue" not in args or args.index("--queue") + 1 >= len(args):
        return False
    q = Path(args[args.index("--queue") + 1])
    q = q if q.is_absolute() else Path(cwd) / q
    try:
        return q.resolve() == Path(queue_db).resolve()
    except OSError:
        return False


def write_outputs(q, root, final=True):
    """Summary files from the queue. Returns (total, unfinished, partial). With
    final=False (a rolling update during the run) only candidates_all.csv and
    done_all.txt are written: the stars still queued aren't failures."""
    rows = q.rows()          # rank, star, state, attempts, worker, status, finished
    final_worker = {star: w for _, star, state, _, w, _, _ in rows if state == "done"}
    # candidates: each star's rows from the worker whose attempt finished it (an
    # earlier failed attempt elsewhere may have left rows in another worker's file)
    files = sorted(root.glob("w*/candidates.csv"), key=lambda p: int(p.parent.name[1:] or 0))
    header, out_rows = None, []
    for f in files:
        slot = int(f.parent.name[1:])
        with open(f, newline="") as fh:
            rd = csv.DictReader(fh)
            if header is None and rd.fieldnames:
                header = rd.fieldnames
            out_rows += [r for r in rd if final_worker.get(r.get("target")) == slot]
    # safety net: one row per (star, candidate), the last written
    out_rows = list({(r.get("target"), r.get("candidate")): r for r in out_rows}.values())
    tmp = root / "candidates_all.csv.tmp"
    with open(tmp, "w", newline="") as fh:
        if header:
            w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore")
            w.writeheader()
            w.writerows(out_rows)
    os.replace(tmp, root / "candidates_all.csv")
    # one line per attempted star: its latest status
    with open(root / "done_all.txt", "w") as fh:
        for _, star, state, att, w, status, finished in rows:
            if status:
                fh.write(f"{star}\t{status}\t{finished or ''}\t{'' if w is None else f'w{w}'}\n")
    unfinished = [star for _, star, state, *_ in rows if state != "done"]
    partial = [(star, status) for _, star, state, _, _, status, _ in rows
               if state == "done" and str(status).startswith("partial_")]
    if not final:
        return len(rows), len(unfinished), len(partial)
    (root / "failed_targets.txt").write_text("".join(f"{s}\n" for s in unfinished))
    (root / "attention.txt").write_text("".join(f"{s}\t{st}\n" for s, st in partial))
    return len(rows), len(unfinished), len(partial)


def rolling_ingest(q, root, py, done, total, state="in progress"):
    """Load the run's results so far into tess.db (tessdb.py ingest replaces the
    run's rows, so repeating it is safe). A failure is logged, not fatal."""
    write_outputs(q, root, final=False)
    note = f"{RUN_NOTE} [{state}: {done}/{total} searched]".strip()
    try:
        r = subprocess.run([py, str(HERE / "tessdb.py"), "ingest", str(root), "--note", note],
                           capture_output=True, text=True, timeout=600, cwd=str(HERE))
        if r.returncode == 0:
            log(f"  rolling update: {done}/{total} stars loaded into the database")
        else:
            log(f"  rolling update failed (exit {r.returncode}): {(r.stderr.strip().splitlines() or ['?'])[-1][:200]}")
    except Exception as e:
        log(f"  rolling update failed: {type(e).__name__}: {e}")


def wait_for_prefetch(root):
    """If prefetch.py is still filling this folder (it holds .prefetch.lock), wait for it,
    so the workers read a complete manifest instead of searching half the stars the
    slow way. Ctrl-C while waiting stops the run as usual."""
    path = root / ".prefetch.lock"
    if not path.exists():
        return
    started, noted = time.time(), None
    with open(path) as fh:
        while True:
            try:
                fcntl.flock(fh, fcntl.LOCK_SH | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if noted is None or time.time() - noted >= 300:
                    log(f"waiting for the prefetch of {root} to finish "
                        f"({(time.time() - started) / 60:.0f} min so far)")
                    noted = time.time()
                time.sleep(PREFETCH_POLL)
        fcntl.flock(fh, fcntl.LOCK_UN)
    if noted is not None:
        log(f"prefetch finished; starting after {(time.time() - started) / 60:.1f} min")


def main():
    if len(sys.argv) < 3:
        print(f"usage: {sys.argv[0]} targets_file out_dir [n_workers] [script]", file=sys.stderr)
        return 1
    targets_file, root = Path(sys.argv[1]), Path(sys.argv[2])
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    script = sys.argv[4] if len(sys.argv) > 4 else "find_planets.py"
    if not targets_file.is_file():
        print(f"no such target file: {targets_file}", file=sys.stderr)
        return 1
    if n < 1:
        print("need at least one worker", file=sys.stderr)
        return 1
    root.mkdir(parents=True, exist_ok=True)
    lock = open(root / ".run.lock", "w")             # one run per output folder
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(f"another run is already using {root}", file=sys.stderr)
        return 4
    stop_file = root / "STOP"
    if stop_file.exists():
        stop_file.unlink()
        log("removed the STOP file left by the last run; continuing")

    targets = read_targets(targets_file)
    (root / "all_targets.txt").write_text("".join(f"{t}\n" for t in targets))
    new_queue = not (root / "queue.db").exists()
    q = Queue(root / "queue.db")
    q.set_meta(max_attempts=MAX_ATTEMPTS, retry_wait=RETRY_WAIT)
    added = q.add(targets)
    if new_queue and list(root.glob("w*/done.txt")):
        log(f"seeded the queue from an earlier run: {seed_from_old_run(q, root)} star(s) already done")
    elif not new_queue and added:
        log(f"{added} new star(s) added to the queue")
    # a rerun gives stars that used up their attempts a fresh set
    fresh = q.con.execute("UPDATE queue SET state = 'todo', attempts = 0, not_before = 0 "
                          "WHERE state = 'failed'").rowcount
    if fresh:
        log(f"{fresh} failed star(s) get fresh attempts")
    commit = write_run_info(root, targets_file, script, n)
    for star, pid, worker in q.running():            # claims of workers from an earlier run
        if not our_worker(pid, root / "queue.db"):
            q.release(star, reason="worker gone (found at start)")
            log(f"  {star}: claim of a dead worker (pid {pid}) returned to the queue")
        else:
            log(f"  {star}: still being searched by pid {pid} from an earlier run")

    py = HERE / ".venv" / "bin" / "python"
    py = str(py) if py.exists() else sys.executable
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               RUN_COMMIT=commit, PYTHONUNBUFFERED="1")
    wait_for_prefetch(root)
    manifest = root / "prefetch.jsonl"           # prefetch.py: stars read from the cache, no MAST
    if manifest.exists():
        env["PREFETCH_MANIFEST"] = str(manifest.resolve())
        n_pre = sum(1 for line in manifest.read_text().splitlines() if '"status": "ok"' in line)
        log(f"prefetch manifest found: {n_pre} star(s) fetched ahead (searched without MAST)")
    c = q.counts()
    log(f"{c['total']} targets, {n} workers, {script} -> {root} (queue: {c['done']} done, "
        f"{c['todo'] + c['running']} to do; up to {MAX_ATTEMPTS} attempts each)")

    interrupted = []
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, lambda s, f: interrupted.append(s))
    children = {}            # slot -> Popen
    last_spawn = -1e9
    last_rolled = c["done"]  # rolling updates count stars finished in this run
    waiting_logged = None
    while True:
        # workers that have exited, and the claims they left behind
        for slot, p in list(children.items()):
            rc = p.poll()
            if rc is None:
                continue
            del children[slot]
            if rc != 0:
                log(f"  worker {slot} (pid {p.pid}) exited with status {rc}: see {root}/w{slot}/stdout.log")
            for star, pid, _ in q.running():
                if pid == p.pid:
                    q.release(star, reason=f"worker exited with status {rc}")
                    log(f"  {star}: returned to the queue (its worker died)")
        ours = {p.pid for p in children.values()}
        for star, pid, _ in q.running():
            if pid not in ours and not our_worker(pid, root / "queue.db"):
                q.release(star, reason="worker gone")
                log(f"  {star}: claim of a dead worker (pid {pid}) returned to the queue")
        if interrupted:
            log(f"interrupted (signal {interrupted[0]}): stopping workers")
            for p in children.values():
                p.terminate()
            for slot, p in children.items():
                try:
                    p.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    p.kill()
                    p.wait()
                for star, pid, _ in q.running():
                    if pid == p.pid:
                        q.release(star, penalise=False)
            children.clear()
            break
        stop = stop_file.exists()
        c = q.counts()
        if ROLLING_INGEST and c["done"] - last_rolled >= ROLLING_INGEST and (c["todo"] or c["running"]):
            rolling_ingest(q, root, py, c["done"], c["total"])
            last_rolled = c["done"]
        if c["todo"] == 0 and c["running"] == 0 and not children:
            break
        if stop and not children:
            break
        free = [k for k in range(n) if k not in children]
        if not stop and c["due"] > 0 and free and time.time() - last_spawn >= STAGGER:
            k = free[0]
            wdir = root / f"w{k}"
            wdir.mkdir(exist_ok=True)
            with open(wdir / "stdout.log", "a") as out:
                p = subprocess.Popen([py, str(HERE / script), "--queue", str(root / "queue.db"),
                                      "--worker", str(k), "--out", str(wdir)],
                                     stdout=out, stderr=subprocess.STDOUT, env=env, cwd=str(HERE))
            children[k] = p
            last_spawn = time.time()
            log(f"  worker {k} started: pid {p.pid} ({c['due']} star(s) due, {c['done']}/{c['total']} done)")
            waiting_logged = None
            continue
        if not children and c["due"] == 0 and c["todo"] > 0 and not stop:
            wait = max(0.0, (q.next_due() or time.time()) - time.time())
            if waiting_logged is None:
                log(f"  {c['todo']} star(s) waiting to be retried; next due in {wait:.0f} s")
                waiting_logged = True
            time.sleep(min(POLL, wait) if wait > 0 else 0)
            continue
        time.sleep(POLL)

    total, nfail, nattn = write_outputs(q, root)
    stopped = bool(interrupted) or stop_file.exists()
    if nfail and stopped:
        log(f"stopped: {total - nfail}/{total} targets done; {nfail} left; needing attention (partial): {nattn}")
    else:
        log(f"complete: {total - nfail}/{total} targets; failed: {nfail}; needing attention (partial): {nattn}")
    if nfail and stopped:
        if ROLLING_INGEST:                       # what's done so far is browsable while paused
            rolling_ingest(q, root, py, total - nfail, total, state="stopped")
        log(f"STOPPED with {nfail} star(s) left; run the same command again to continue")
        return 6
    if nfail:
        log(f"INCOMPLETE: unfinished stars in {root}/failed_targets.txt; run the same command again to retry")
        return 3
    for f in ("failed_targets.txt", "attention.txt"):
        if (root / f).exists() and not (root / f).stat().st_size:
            (root / f).unlink()
    return 0


if __name__ == "__main__":
    sys.exit(main())
