#!/usr/bin/env python3
"""Timestamped fingerprint of a candidate, without revealing it.

Writes investigations/tic<id>/discovery_record.txt (star, signals, batch, date found),
in a fixed plain-text format, and appends a line to the public PlanetCandidates
discovery_log.md: the date found, a neutral description, and the SHA-256 of the record.
The record stays in the private repo. Revealing it later (for example when the
candidate is posted) lets anyone check that its hash matches the dated public line.

Usage: .venv/bin/python investigations/fingerprint.py "TIC 123456789" RUN "description" [CAND ...]
    RUN          the batch the candidate was found in (e.g. b022)
    description  neutral, e.g. "three transit signals on one K dwarf"
    CAND         candidate numbers to include (default: every non-same_period signal)
The log is ~/candidates/discovery_log.md (PLANETCANDIDATES_LOG to override).
FINGERPRINT_SOURCE=search uses the search's ephemeris instead of the transit fit's. Never
edit a record after its hash is logged: the hash would no longer match.
"""
import hashlib
import json
import os
import re
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = Path(os.environ.get("PLANETCANDIDATES_LOG", Path.home() / "candidates" / "discovery_log.md"))
HEADER = """# Discovery log

Dated fingerprints of candidates found by PlanetFinder, an independent search of TESS
light curves, logged before they are made public. Each line is the SHA-256 hash of a
plain-text discovery record (star, signals, batch, date) kept privately. Revealing a
record later lets anyone check that its hash matches the line below, which shows the
candidate was found by that date. This file's git history on GitHub dates each line.

| Logged | Found (batch) | Description | SHA-256 of discovery record | Record |
|---|---|---|---|---|
"""


def record(con, star, run, cands):
    con.row_factory = sqlite3.Row
    r = con.execute("SELECT run_id, started FROM runs WHERE name = ?", (run,)).fetchone()
    if not r:
        sys.exit(f"no run {run}")
    s = con.execute("SELECT * FROM stars WHERE star = ?", (star,)).fetchone()
    rows = con.execute("SELECT * FROM candidates WHERE run_id = ? AND star = ? AND status != 'same_period' "
                       "ORDER BY candidate", (r["run_id"], star)).fetchall()
    rows = [x for x in rows if not cands or x["candidate"] in cands]
    if not rows:
        sys.exit(f"no signals for {star} in {run}")
    lines = ["PlanetFinder discovery record",
             f"star: {star}",
             f"position_deg: {s['ra']:.5f} {s['dec']:+.5f}" if s and s["ra"] is not None else "position_deg: unknown",
             f"tmag: {s['tmag']:.2f}" if s and s["tmag"] is not None else "tmag: unknown",
             f"found_in: batch {run}, search started {r['started'][:16]}"]
    for x in rows:
        fit = con.execute("SELECT data FROM transit_fits WHERE run_id = ? AND star = ? AND candidate = ?",
                          (r["run_id"], star, x["candidate"])).fetchone()
        if fit and os.environ.get("FINGERPRINT_SOURCE") != "search":   # search: when the fit went wrong
            f = json.loads(fit[0])
            P, T0, dep, dur = (float(f[k]) for k in ("period", "t0", "depth_ppm", "t14_h"))
            src = "transit fit"
        else:
            P, T0, dep, dur, src = x["period_d"], x["t0_btjd"], x["depth_ppm"], x["duration_h"], "search"
        lines.append(f"signal {x['candidate']}: period_d={P:.6f} t0_btjd={T0:.5f} depth_ppm={dep:.0f} "
                     f"duration_h={dur:.2f} rp_earth={x['rp_earth']} snr={x['snr']} ({src})")
    lines.append("")
    return "\n".join(lines)


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    star, run, desc = sys.argv[1], sys.argv[2], sys.argv[3]
    cands = [int(c) for c in sys.argv[4:]]
    tic = re.sub(r"\D", "", star)
    path = ROOT / "investigations" / f"tic{tic}" / "discovery_record.txt"
    if path.exists():
        sys.exit(f"{path} already exists (and may already be logged): not overwritten")
    con = sqlite3.connect(f"file:{ROOT / 'db' / 'tess.db'}?mode=ro", uri=True)
    text = record(con, star, run, cands)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    digest = hashlib.sha256(text.encode()).hexdigest()
    found = con.execute("SELECT substr(started, 1, 10) FROM runs WHERE name = ?", (run,)).fetchone()[0]
    import datetime
    if not LOG.exists():
        LOG.write_text(HEADER)
    with open(LOG, "a") as fh:
        fh.write(f"| {datetime.date.today()} | {found} ({run}) | {desc} | `{digest}` | not yet revealed |\n")
    print(f"{path}\nsha256 {digest}\nlogged in {LOG}")


if __name__ == "__main__":
    main()
