"""
workqueue.py - the shared target queue of a run: <out>/queue.db (SQLite).

One row per star, in target-file order (rank). Workers claim the best-ranked
star that is due, search it and record the outcome; a star that ends in an error
goes back into the queue after a delay, until it has used its attempts. The
supervisor (run_parallel.py) creates the queue, returns the claims of workers
that died, and writes the run's summary files from it.

States: todo -> running -> done | todo (an error, retried later) | failed
(an error on the last attempt). Claiming is one transaction (BEGIN IMMEDIATE),
so two workers can never get the same star.
"""
import datetime
import os
import sqlite3
import time

SCHEMA = """
CREATE TABLE IF NOT EXISTS queue (
    rank INTEGER PRIMARY KEY,          -- position in the target file (score order)
    star TEXT UNIQUE NOT NULL,
    state TEXT NOT NULL DEFAULT 'todo',
    attempts INTEGER NOT NULL DEFAULT 0,
    not_before REAL NOT NULL DEFAULT 0, -- unix time a retry is due
    worker INTEGER, pid INTEGER, commit_sha TEXT,
    claimed TEXT, finished TEXT, status TEXT);
CREATE INDEX IF NOT EXISTS queue_state ON queue (state, not_before, rank);
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
"""


def now_iso():
    return datetime.datetime.now().isoformat(timespec="seconds")


def is_error(status):
    return str(status).startswith("error")


class Queue:
    def __init__(self, path):
        self.path = str(path)
        self.con = sqlite3.connect(self.path, timeout=60, isolation_level=None)
        self.con.execute("PRAGMA journal_mode=WAL")
        self.con.execute("PRAGMA busy_timeout=60000")
        self.con.executescript(SCHEMA)

    # ---------------------------------------------------------------- settings
    def set_meta(self, **kv):
        for k, v in kv.items():
            self.con.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)", (k, str(v)))

    def meta(self, key, default=None):
        r = self.con.execute("SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
        return r[0] if r else default

    # ---------------------------------------------------------------- filling
    def add(self, stars):
        """Add stars not yet queued, after the existing ones, in the given order."""
        have = {r[0] for r in self.con.execute("SELECT star FROM queue")}
        top = self.con.execute("SELECT coalesce(max(rank), 0) FROM queue").fetchone()[0]
        new = [s for s in dict.fromkeys(stars) if s not in have]
        self.con.execute("BEGIN IMMEDIATE")
        self.con.executemany("INSERT INTO queue (rank, star) VALUES (?, ?)",
                             [(top + i + 1, s) for i, s in enumerate(new)])
        self.con.execute("COMMIT")
        return len(new)

    def mark_done(self, star, status, finished=None, worker=None):
        """Record an outcome from outside a claim (seeding from an old run's done.txt)."""
        self.con.execute("UPDATE queue SET state = 'done', status = ?, finished = ?, worker = ? "
                         "WHERE star = ? AND state != 'done'", (status, finished or now_iso(), worker, star))

    # ---------------------------------------------------------------- workers
    def claim(self, worker, pid=None, commit=None):
        """The best-ranked star that is due, now marked running for this worker;
        None when nothing is due."""
        self.con.execute("BEGIN IMMEDIATE")
        try:
            r = self.con.execute("SELECT rank, star FROM queue WHERE state = 'todo' AND not_before <= ? "
                                 "ORDER BY rank LIMIT 1", (time.time(),)).fetchone()
            if r:
                self.con.execute("UPDATE queue SET state = 'running', attempts = attempts + 1, worker = ?, "
                                 "pid = ?, commit_sha = ?, claimed = ? WHERE rank = ?",
                                 (worker, pid or os.getpid(), commit, now_iso(), r[0]))
            self.con.execute("COMMIT")
        except BaseException:
            self.con.execute("ROLLBACK")
            raise
        return r[1] if r else None

    def finish(self, star, status):
        """Record a claimed star's outcome. An error is retried after retry_wait
        seconds, unless that was its last attempt."""
        max_att = int(self.meta("max_attempts", 3))
        wait = float(self.meta("retry_wait", 300))
        self.con.execute("BEGIN IMMEDIATE")
        att = self.con.execute("SELECT attempts FROM queue WHERE star = ?", (star,)).fetchone()[0]
        if not is_error(status):
            state, nb = "done", 0
        elif att >= max_att:
            state, nb = "failed", 0
        else:
            state, nb = "todo", time.time() + wait
        self.con.execute("UPDATE queue SET state = ?, status = ?, finished = ?, not_before = ?, pid = NULL "
                         "WHERE star = ?", (state, status, now_iso(), nb, star))
        self.con.execute("COMMIT")
        return state

    # ---------------------------------------------------------------- supervisor
    def release(self, star, penalise=True, reason="worker died"):
        """Return a running star to the queue (its worker is gone). penalise=False
        undoes the attempt (the run was stopped, the star wasn't at fault)."""
        max_att = int(self.meta("max_attempts", 3))
        wait = float(self.meta("retry_wait", 300))
        self.con.execute("BEGIN IMMEDIATE")
        r = self.con.execute("SELECT attempts FROM queue WHERE star = ? AND state = 'running'", (star,)).fetchone()
        if r:
            att = r[0] if penalise else max(0, r[0] - 1)
            if penalise and att >= max_att:
                self.con.execute("UPDATE queue SET state = 'failed', attempts = ?, pid = NULL, status = ?, "
                                 "finished = ? WHERE star = ?", (att, f"error: {reason}", now_iso(), star))
            else:
                note = f"error: {reason}" if penalise else None
                self.con.execute("UPDATE queue SET state = 'todo', attempts = ?, pid = NULL, not_before = ?, "
                                 "status = CASE WHEN ? IS NULL THEN status ELSE ? END WHERE star = ?",
                                 (att, time.time() + (wait if penalise else 0), note, note, star))
        self.con.execute("COMMIT")
        return bool(r)

    def running(self):
        return self.con.execute("SELECT star, pid, worker FROM queue WHERE state = 'running'").fetchall()

    def counts(self):
        c = {"todo": 0, "running": 0, "done": 0, "failed": 0, "due": 0}
        for state, n in self.con.execute("SELECT state, count(*) FROM queue GROUP BY state"):
            c[state] = n
        c["due"] = self.con.execute("SELECT count(*) FROM queue WHERE state = 'todo' AND not_before <= ?",
                                    (time.time(),)).fetchone()[0]
        c["total"] = sum(c[s] for s in ("todo", "running", "done", "failed"))
        return c

    def next_due(self):
        r = self.con.execute("SELECT min(not_before) FROM queue WHERE state = 'todo'").fetchone()[0]
        return r

    def rows(self):
        return self.con.execute("SELECT rank, star, state, attempts, worker, status, finished "
                                "FROM queue ORDER BY rank").fetchall()

    def close(self):
        self.con.close()
