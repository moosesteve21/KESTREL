#!/usr/bin/env python3
"""
tessdb.py - central SQLite record of PlanetFinder runs, candidates and notes.

Usage:
    python tessdb.py ingest results/aries [--name aries] [--commit SHA] [--note TEXT]
    python tessdb.py note "TIC 434331752" "VSX submitted; provisional ATO J093.6734+17.5725" --tag vsx
    python tessdb.py list                      # runs in the database
    python tessdb.py remove RUN                # delete a run (e.g. a test)
    python tessdb.py backup                    # dated copy in db/backups/
    python tessdb.py clean-cache [--dry-run]   # delete cached light curves of no-signal stars
    python tessdb.py serve [--port 8001]       # browse with Datasette on 127.0.0.1

Runs are ingested after they finish; the per-worker CSVs stay as the raw record.
Re-ingesting a run replaces its rows, so it is safe to repeat. A folder can hold
a search run (candidates_all.csv or candidates.csv), enriched.csv, pixel-check
output (pixel_check*/), reduction-check output (reduction_check*/) and/or
injections.csv. run_info.json (written by
run_parallel.sh) gives the code version and settings; for older folders pass
--commit and --note.

The database file is db/tess.db (TESSDB environment variable to override);
dated backups go to db/backups/ and Datasette's settings to db/datasette/. It is not in git: back it up with `python tessdb.py backup`.
"""
import argparse
import csv
import datetime
import json
import os
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DB = Path(os.environ.get("TESSDB", HERE / "db" / "tess.db"))
BACKUPS = HERE / "db" / "backups" if "TESSDB" not in os.environ else DB.parent / "backups"   # scratch DBs back up beside themselves
DATASETTE = HERE / "db" / "datasette"
# records/: the tables nothing else can recreate, as CSV in git (written by `backup`).
# A scratch database (TESSDB set) keeps its records beside it, out of the repo.
RECORDS = HERE / "records" if "TESSDB" not in os.environ else DB.parent / "records"
RECORD_TABLES = {                  # table -> sort order (stable, so git diffs stay small)
    "notes": "created, star, tag, text",
    "watchlist": "star",
    "watch_sectors": "star, sector",
    "watch_events": "star, sector, epoch",
    "batches": "name",
}
PORT = 8001

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    run_id INTEGER PRIMARY KEY, name TEXT UNIQUE, folder TEXT, kind TEXT,
    started TEXT, ingested TEXT, commit_sha TEXT, uncommitted_changes INTEGER,
    targets_file TEXT, n_targets INTEGER, settings TEXT, note TEXT);
CREATE TABLE IF NOT EXISTS stars (
    star TEXT PRIMARY KEY, tic TEXT, ra REAL, dec REAL, tmag REAL, teff REAL,
    rad REAL, mass REAL, contratio REAL, n_sectors INTEGER, sectors TEXT, field TEXT);
CREATE TABLE IF NOT EXISTS searches (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, star TEXT, status TEXT,
    finished TEXT, worker TEXT, PRIMARY KEY (run_id, star));
CREATE TABLE IF NOT EXISTS candidates (
    cand_id INTEGER PRIMARY KEY, run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE,
    star TEXT, candidate INTEGER, period_d REAL, t0_btjd REAL, duration_h REAL,
    depth_ppm REAL, snr REAL, sde REAL, dips TEXT, n_transits INTEGER, rp_earth REAL,
    status TEXT, flags TEXT, enrich_flags TEXT, plot TEXT, data TEXT,
    UNIQUE (run_id, star, candidate));
CREATE TABLE IF NOT EXISTS candidate_flags (
    cand_id INTEGER REFERENCES candidates(cand_id) ON DELETE CASCADE, flag TEXT, source TEXT);
CREATE TABLE IF NOT EXISTS pixel_checks (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, star TEXT, candidate INTEGER,
    sector INTEGER, depth_1px_ppm REAL, depth_3x3_ppm REAL, peak_sigma REAL,
    offset_arcsec REAL, data TEXT, folder TEXT);
CREATE TABLE IF NOT EXISTS pixel_verdicts (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, star TEXT, candidate INTEGER,
    verdict TEXT, source TEXT, folder TEXT);
CREATE TABLE IF NOT EXISTS depth_maps (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, star TEXT, candidate INTEGER,
    verdict TEXT, source TEXT, target_dchi2 REAL, detect_dchi2 REAL, sectors_won TEXT,
    ephemeris TEXT, detail TEXT, folder TEXT);
CREATE TABLE IF NOT EXISTS reduction_checks (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, star TEXT, candidate INTEGER,
    verdict TEXT, n_support INTEGER, n_contradict INTEGER, reference_ppm REAL, detail TEXT,
    folder TEXT);
CREATE TABLE IF NOT EXISTS transit_fits (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, star TEXT, candidate INTEGER,
    verdict TEXT, density_z REAL, rho REAL, rho_lo REAL, rho_hi REAL, rho_tic REAL,
    rp_earth REAL, rp_rs REAL, b REAL, a_rs REAL, t14_h REAL, period REAL, t0 REAL,
    n_events INTEGER, data TEXT, folder TEXT);
CREATE TABLE IF NOT EXISTS fpp_checks (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, star TEXT, candidate INTEGER,
    label TEXT, fpp REAL, nfpp REAL, prob_tp REAL, draws INTEGER, sectors TEXT,
    scenarios TEXT, detail TEXT, folder TEXT);
CREATE TABLE IF NOT EXISTS injections (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, job INTEGER, star TEXT,
    p_inj REAL, rp_inj REAL, depth_inj_ppm REAL, snr_exp REAL, recovered TEXT,
    rec_status TEXT, rec_flags TEXT, data TEXT);
CREATE TABLE IF NOT EXISTS vsx_shortlist (
    run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, star TEXT, candidate INTEGER,
    grp TEXT, signature TEXT, vsx TEXT, other_catalogues TEXT, simbad_type TEXT,
    blend_suspects TEXT, UNIQUE (run_id, star, candidate));
CREATE TABLE IF NOT EXISTS watchlist (
    star TEXT PRIMARY KEY, period REAL, period_err REAL, t0 REAL, t0_err REAL,
    duration_h REAL, depth_ppm REAL, baseline_sectors TEXT, reason TEXT, source TEXT,
    added TEXT, last_checked TEXT, verdict TEXT, summary TEXT, next_sectors TEXT);
CREATE TABLE IF NOT EXISTS watch_sectors (
    star TEXT REFERENCES watchlist(star) ON DELETE CASCADE, sector INTEGER, product TEXT,
    baseline INTEGER, n_events INTEGER, checked TEXT, PRIMARY KEY (star, sector));
CREATE TABLE IF NOT EXISTS watch_events (
    star TEXT REFERENCES watchlist(star) ON DELETE CASCADE, sector INTEGER, epoch INTEGER,
    tc REAL, depth_ppm REAL, err_ppm REAL, snr REAL, baseline INTEGER, product TEXT,
    checked TEXT, PRIMARY KEY (star, sector, epoch));
CREATE TABLE IF NOT EXISTS notes (
    note_id INTEGER PRIMARY KEY, created TEXT, star TEXT, run TEXT, candidate INTEGER,
    tag TEXT, text TEXT);
CREATE TABLE IF NOT EXISTS catalog (
    tic TEXT PRIMARY KEY, ra REAL, dec REAL, tmag REAL, teff REAL, rad REAL, mass REAL,
    logg REAL, contratio REAL, gaia TEXT, glat REAL, elat REAL, n_sectors INTEGER,
    sectors TEXT, known TEXT, ruwe REAL, score REAL, batch TEXT, updated TEXT);
CREATE TABLE IF NOT EXISTS catalog_tiles (tile TEXT PRIMARY KEY, n INTEGER, done TEXT);
CREATE TABLE IF NOT EXISTS batches (
    name TEXT PRIMARY KEY, created TEXT, n INTEGER, criteria TEXT, targets_file TEXT);
CREATE INDEX IF NOT EXISTS i_cat_score ON catalog (score);
CREATE INDEX IF NOT EXISTS i_cand_star ON candidates (star);
CREATE INDEX IF NOT EXISTS i_flags ON candidate_flags (flag);
-- one row per check: a duplicate in a run's output is an error, not something to keep
CREATE UNIQUE INDEX IF NOT EXISTS u_pixel_checks ON pixel_checks (run_id, star, candidate, sector);
CREATE UNIQUE INDEX IF NOT EXISTS u_pixel_verdicts ON pixel_verdicts (run_id, star, candidate);
CREATE UNIQUE INDEX IF NOT EXISTS u_reduction_checks ON reduction_checks (run_id, star, candidate);
CREATE UNIQUE INDEX IF NOT EXISTS u_transit_fits ON transit_fits (run_id, star, candidate);
CREATE UNIQUE INDEX IF NOT EXISTS u_fpp_checks ON fpp_checks (run_id, star, candidate);
CREATE UNIQUE INDEX IF NOT EXISTS u_depth_maps ON depth_maps (run_id, star, candidate);
CREATE UNIQUE INDEX IF NOT EXISTS u_injections ON injections (run_id, job);
DROP VIEW IF EXISTS candidates_view;
CREATE VIEW candidates_view AS
    SELECT c.cand_id, r.name AS run, c.star, c.candidate, c.status,
           round(c.period_d, 4) AS period_d, c.depth_ppm, c.snr, c.rp_earth,
           CASE WHEN c.plot IS NULL OR c.plot = '' THEN NULL
                WHEN r.folder LIKE '%/results/' || r.name
                THEN 'http://localhost:{port}/results/' || r.name || '/' || c.plot
                ELSE 'http://localhost:{port}/' || r.name || '/' || c.plot END AS plot_url,
           CASE WHEN instr(c.plot, '/') > 0 AND substr(c.plot, 1, 1) = 'w'
                THEN substr(c.plot, 1, instr(c.plot, '/') - 1) END AS worker,
           c.flags, c.enrich_flags,
           (SELECT f.verdict || ' (' || printf('%+.1f', f.density_z) || ' sigma)' FROM transit_fits f
            WHERE f.run_id = c.run_id AND f.star = c.star AND f.candidate = c.candidate) AS density,
           (SELECT d.verdict || CASE WHEN d.verdict = 'off_target' THEN ': ' || d.source ELSE '' END
            FROM depth_maps d
            WHERE d.run_id = c.run_id AND d.star = c.star AND d.candidate = c.candidate) AS depth_map,
           (SELECT x.label || CASE WHEN x.fpp IS NULL THEN '' ELSE ' (FPP ' || printf('%.3f', x.fpp) || ')' END
            FROM fpp_checks x
            WHERE x.run_id = c.run_id AND x.star = c.star AND x.candidate = c.candidate) AS fpp,
           (SELECT group_concat(n.tag || ': ' || n.text, ' | ') FROM notes n WHERE n.star = c.star)
               AS notes
    FROM candidates c JOIN runs r USING (run_id);
DROP VIEW IF EXISTS vsx_shortlist_view;
CREATE VIEW vsx_shortlist_view AS
    SELECT v.grp AS "group", r.name AS run, v.star, v.candidate, c.period_d, c.depth_ppm,
           c.snr, v.signature, v.vsx, v.other_catalogues, v.simbad_type, v.blend_suspects,
           c.plot_url, c.notes, c.flags, c.enrich_flags
    FROM vsx_shortlist v JOIN runs r USING (run_id)
         LEFT JOIN candidates_view c ON c.run = r.name AND c.star = v.star AND c.candidate = v.candidate;
DROP VIEW IF EXISTS stars_searched;
CREATE VIEW stars_searched AS
    SELECT s.star, count(*) AS runs, group_concat(r.name, ', ') AS run_names,
           max(s.finished) AS last_searched,
           (SELECT count(*) FROM candidates c WHERE c.star = s.star) AS signals
    FROM searches s JOIN runs r USING (run_id) GROUP BY s.star;
"""

SEARCH_COLS = ("period_d", "t0_btjd", "duration_h", "depth_ppm", "snr", "sde", "dips",
               "n_transits", "rp_earth", "status", "flags", "plot")


def connect():
    con = sqlite3.connect(DB, timeout=120)   # a search's rolling ingest and another batch's follow-up may overlap
    con.execute("PRAGMA foreign_keys = ON")
    con.executescript(SCHEMA.replace("{port}", str(PORT)))
    return con


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def resolve_plot(folder, plot):
    """Plot path relative to the run folder. Workers write paths relative to
    their own folder, and the merged CSV doesn't say which worker: find it."""
    if not plot:
        return ""
    if (folder / plot).exists():
        return plot
    hits = sorted(folder.glob(f"w*/{plot}"))
    return str(hits[0].relative_to(folder)) if hits else plot


def find_targets_csv(targets_file):
    """The finder's CSV for a target list (targets/aries.txt -> targets/aries.csv),
    for star properties."""
    if not targets_file:
        return None
    p = (HERE / targets_file) if not Path(targets_file).is_absolute() else Path(targets_file)
    cand = p.with_suffix(".csv")
    return cand if cand.exists() else None


def add(con, where, sql, values):
    """INSERT one row; a row that breaks a uniqueness rule stops the ingest (the
    transaction is rolled back, so the database keeps its previous state)."""
    try:
        return con.execute(sql, values)
    except sqlite3.IntegrityError as e:
        sys.exit(f"{e}; each check may appear once per run: remove the duplicate (or the extra folder)\n"
                 f"ingest refused, database unchanged: duplicate row in {where}")


def load_vsx_shortlist(con, run_id, folder):
    """The run's vsx_shortlist.csv (vsx_shortlist.py), if it has one."""
    f = Path(folder) / "vsx_shortlist.csv"
    if not f.exists():
        return
    for r in read_csv(f):
        add(con, f"vsx_shortlist.csv: {r['star']} candidate {r['candidate']}",
            "INSERT INTO vsx_shortlist VALUES (?,?,?,?,?,?,?,?,?)",
            (run_id, r["star"], int(r["candidate"]), r["group"], r["signature"], r["vsx"],
             r["other_catalogues"], r["simbad_type"], r["blend_suspects"]))


def ingest(folder, name=None, commit=None, note=None, targets=None, started=None):
    folder = Path(folder).resolve()
    if not folder.is_dir():
        sys.exit(f"no such folder: {folder}")
    name = name or folder.name
    info = {}
    if (folder / "run_info.json").exists():
        info = json.loads((folder / "run_info.json").read_text())
    if targets:
        info["targets_file"] = targets
    search_csv = next((folder / f for f in ("candidates_all.csv", "candidates.csv")
                       if (folder / f).exists()), None)
    done_files = ([folder / "done_all.txt"] if (folder / "done_all.txt").exists()
                  else sorted(folder.glob("w*/done.txt")) or
                  ([folder / "done.txt"] if (folder / "done.txt").exists() else []))
    kinds = [k for k, ok in (("search", search_csv), ("injection", (folder / "injections.csv").exists()))
             if ok]
    if not kinds:
        sys.exit(f"{folder}: no candidates CSV or injections.csv found")
    started = started or info.get("started")
    if not started:                                     # earliest log line, else folder time
        logs = sorted(folder.glob("w*/run.log")) + sorted(folder.glob("run.log"))
        stamps = [l.read_text().split(" INFO", 1)[0][:19] for l in logs if l.stat().st_size]
        started = min(stamps) if stamps else datetime.datetime.fromtimestamp(
            folder.stat().st_mtime).isoformat(timespec="seconds")
    con = connect()
    with con:
        con.execute("DELETE FROM runs WHERE name = ?", (name,))      # replace on re-ingest
        cur = con.execute(
            "INSERT INTO runs (name, folder, kind, started, ingested, commit_sha, uncommitted_changes,"
            " targets_file, settings, note) VALUES (?,?,?,?,?,?,?,?,?,?)",
            (name, str(folder), "+".join(kinds), started.replace("T", " "),
             datetime.datetime.now().isoformat(timespec="seconds"),
             commit or info.get("commit"), info.get("uncommitted_changes"),
             info.get("targets_file"), json.dumps(info.get("settings")) if info.get("settings") else None,
             note))
        run_id = cur.lastrowid
        # stars: from the finder's CSV when the target list came from it
        tcsv = find_targets_csv(info.get("targets_file"))
        field = Path(info.get("targets_file", "")).stem or None
        if tcsv:
            for r in read_csv(tcsv):
                star = r.get("target") or f"TIC {r['ID']}"
                con.execute("INSERT INTO stars (star, tic, ra, dec, tmag, teff, rad, mass, contratio,"
                            " n_sectors, sectors, field) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
                            " ON CONFLICT (star) DO UPDATE SET tmag=excluded.tmag, sectors=excluded.sectors",
                            (star, str(r["ID"]), num(r["ra"]), num(r["dec"]), num(r["Tmag"]),
                             num(r["Teff"]), num(r["rad"]), num(r["mass"]), num(r["contratio"]),
                             num(r["n_sectors"]), r["sectors"], field))
        n_targets = 0
        for f in done_files:
            for line in f.read_text().splitlines():
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                n_targets += 1
                con.execute("INSERT OR REPLACE INTO searches VALUES (?,?,?,?,?)",
                            (run_id, parts[0], parts[1], parts[2] if len(parts) > 2 else None,
                             parts[3] if len(parts) > 3 else       # queue runs: done_all.txt has it
                             f.parent.name if f.parent != folder else ""))
                m = re.fullmatch(r"TIC\s*(\d+)", parts[0])
                con.execute("INSERT OR IGNORE INTO stars (star, tic) VALUES (?,?)",
                            (parts[0], m.group(1) if m else None))
        con.execute("UPDATE runs SET n_targets = ? WHERE run_id = ?", (n_targets, run_id))
        # candidates (+ enrichment where present)
        if search_csv:
            enriched = {}
            if (folder / "enriched.csv").exists():
                enriched = {(r["target"], r["candidate"]): r for r in read_csv(folder / "enriched.csv")}
            for r in read_csv(search_csv):
                e = enriched.get((r["target"], r["candidate"]), {})
                extra = {k: v for k, v in {**r, **e}.items()
                         if k not in SEARCH_COLS and k not in ("target", "candidate") and v not in ("", None)}
                cur = add(con, f"{search_csv.name}: {r['target']} candidate {r['candidate']}",
                    "INSERT INTO candidates (run_id, star, candidate, period_d, t0_btjd, duration_h,"
                    " depth_ppm, snr, sde, dips, n_transits, rp_earth, status, flags, enrich_flags,"
                    " plot, data) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (run_id, r["target"], int(r["candidate"]), num(r["period_d"]), num(r["t0_btjd"]),
                     num(r["duration_h"]), num(r["depth_ppm"]), num(r["snr"]), num(r["sde"]),
                     r.get("dips"), num(r.get("n_transits")), num(r.get("rp_earth")), r["status"],
                     r.get("flags", ""), e.get("enrich_flags", ""), resolve_plot(folder, r.get("plot", "")),
                     json.dumps(extra)))
                cid = cur.lastrowid
                for src, field_ in (("search", r.get("flags", "")), ("enrichment", e.get("enrich_flags", ""))):
                    for fl in filter(None, (field_ or "").split(";")):
                        con.execute("INSERT INTO candidate_flags VALUES (?,?,?)", (cid, fl, src))
        # pixel checks: any pixel_check*/ folder in the run
        for pdir in sorted(folder.glob("pixel_check*")):
            rel = pdir.name
            if (pdir / "pixel_check.csv").exists():
                for r in read_csv(pdir / "pixel_check.csv"):
                    add(con, f"{rel}/pixel_check.csv: {r['target']} candidate {r['candidate']} sector {r['sector']}",
                        "INSERT INTO pixel_checks VALUES (?,?,?,?,?,?,?,?,?,?)",
                                (run_id, r["target"], int(r["candidate"]), int(r["sector"]),
                                 num(r["depth_1px_ppm"]), num(r["depth_3x3_ppm"]),
                                 num(r["peak_sigma"]), num(r["offset_arcsec"]), json.dumps(r), rel))
            if (pdir / "pixel_summary.csv").exists():
                for r in read_csv(pdir / "pixel_summary.csv"):
                    add(con, f"{rel}/pixel_summary.csv: {r['target']} candidate {r['candidate']}",
                        "INSERT INTO pixel_verdicts VALUES (?,?,?,?,?,?)",
                                (run_id, r["target"], int(r["candidate"]), r["verdict"],
                                 r.get("source", ""), rel))
        # reduction checks: any reduction_check*/ folder in the run
        for rdir in sorted(folder.glob("reduction_check*")):
            if not (rdir / "reduction_summary.csv").exists():
                continue
            for r in read_csv(rdir / "reduction_summary.csv"):
                add(con, f"{rdir.name}/reduction_summary.csv: {r['target']} candidate {r['candidate']}",
                    "INSERT INTO reduction_checks VALUES (?,?,?,?,?,?,?,?,?)",
                            (run_id, r["target"], int(r["candidate"]), r["verdict"],
                             num(r.get("n_support")), num(r.get("n_contradict")),
                             num(r.get("reference_ppm")), r.get("detail", ""), rdir.name))
                if r["verdict"] == "reduction_dependent":
                    row = con.execute("SELECT cand_id FROM candidates WHERE run_id = ? AND star = ? "
                                      "AND candidate = ?", (run_id, r["target"], int(r["candidate"]))).fetchone()
                    if row:
                        con.execute("INSERT INTO candidate_flags VALUES (?,?,?)",
                                    (row[0], "reduction_dependent", "reduction_check"))
        # physical transit fits (transit_fit.py) and false-positive probabilities (fpp_check.py)
        for fdir in sorted(folder.glob("transit_fit*")):
            if not (fdir / "transit_fit.csv").exists():
                continue
            for r in read_csv(fdir / "transit_fit.csv"):
                add(con, f"{fdir.name}/transit_fit.csv: {r['target']} candidate {r['candidate']}",
                    "INSERT INTO transit_fits VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            (run_id, r["target"], int(r["candidate"]), r["verdict"], num(r.get("density_z")),
                             num(r.get("rho")), num(r.get("rho_lo")), num(r.get("rho_hi")), num(r.get("rho_tic")),
                             num(r.get("rp_earth")), num(r.get("rp_rs")), num(r.get("b")), num(r.get("a_rs")),
                             num(r.get("t14_h")), num(r.get("period")), num(r.get("t0")),
                             num(r.get("n_events")), json.dumps(r), fdir.name))
                if r["verdict"] in ("density_mismatch", "duration_too_long"):
                    row = con.execute("SELECT cand_id FROM candidates WHERE run_id = ? AND star = ? "
                                      "AND candidate = ?", (run_id, r["target"], int(r["candidate"]))).fetchone()
                    if row:
                        con.execute("INSERT INTO candidate_flags VALUES (?,?,?)",
                                    (row[0], f"fit_{r['verdict']}" if r["verdict"] == "duration_too_long"
                                     else r["verdict"], "transit_fit"))
        for ddir in sorted(folder.glob("depth_map*")):
            if not (ddir / "depth_map.csv").exists():
                continue
            for r in read_csv(ddir / "depth_map.csv"):
                add(con, f"{ddir.name}/depth_map.csv: {r['target']} candidate {r['candidate']}",
                    "INSERT INTO depth_maps VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    (run_id, r["target"], int(r["candidate"]), r["verdict"], r.get("source", ""),
                     num(r.get("target_dchi2")), num(r.get("detect_dchi2")), r.get("sectors_won", ""),
                     r.get("ephemeris", ""), r.get("detail", ""), ddir.name))
                if r["verdict"] == "off_target":
                    row = con.execute("SELECT cand_id FROM candidates WHERE run_id = ? AND star = ? "
                                      "AND candidate = ?", (run_id, r["target"], int(r["candidate"]))).fetchone()
                    if row:
                        con.execute("INSERT INTO candidate_flags VALUES (?,?,?)",
                                    (row[0], "depth_map_off_target", "depth_map"))
        for xdir in sorted(folder.glob("fpp_check*")):
            if not (xdir / "fpp.csv").exists():
                continue
            for r in read_csv(xdir / "fpp.csv"):
                add(con, f"{xdir.name}/fpp.csv: {r['target']} candidate {r['candidate']}",
                    "INSERT INTO fpp_checks VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                            (run_id, r["target"], int(r["candidate"]), r.get("label"), num(r.get("fpp")),
                             num(r.get("nfpp")), num(r.get("prob_tp")), num(r.get("draws")),
                             r.get("sectors", ""), r.get("scenarios", ""), r.get("detail", ""), xdir.name))
        load_vsx_shortlist(con, run_id, folder)
        # injections
        if (folder / "injections.csv").exists():
            for r in read_csv(folder / "injections.csv"):
                add(con, f"injections.csv: job {r['job']} ({r['star']})",
                    "INSERT INTO injections VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                            (run_id, int(r["job"]), r["star"], num(r.get("p_inj")), num(r.get("rp_inj")),
                             num(r.get("depth_inj_ppm")), num(r.get("snr_exp")), r.get("recovered"),
                             r.get("rec_status"), r.get("rec_flags"), json.dumps(r)))
    counts = {t: con.execute(f"SELECT count(*) FROM {t} WHERE run_id = ?", (run_id,)).fetchone()[0]
              for t in ("searches", "candidates", "pixel_checks", "pixel_verdicts", "reduction_checks",
                        "transit_fits", "depth_maps", "fpp_checks", "vsx_shortlist", "injections")}
    print(f"ingested {name} ({'+'.join(kinds)}): " + ", ".join(f"{v} {k}" for k, v in counts.items() if v))
    if not (commit or info.get("commit")):
        print("  note: no code version recorded (no run_info.json); pass --commit if known")


def rebuild():
    """Recreate the database with the current schema: re-ingest every run from its
    folder with its recorded details, and restore the notes. Use after a schema
    change. The old file is kept in db/backups/."""
    global DB
    old = sqlite3.connect(DB)
    old.row_factory = sqlite3.Row
    runs = [dict(r) for r in old.execute("SELECT * FROM runs ORDER BY started")]
    notes = [dict(r) for r in old.execute("SELECT * FROM notes ORDER BY note_id")]
    kept = {}                                        # not rebuilt from run folders: copy across
    for t in ("catalog", "catalog_tiles", "batches", "watchlist", "watch_sectors", "watch_events"):
        try:
            kept[t] = [dict(r) for r in old.execute(f"SELECT * FROM {t}")]
        except sqlite3.OperationalError:
            kept[t] = []
    old.close()
    final, DB = DB, DB.with_name(DB.name + ".new")
    DB.unlink(missing_ok=True)
    for r in runs:
        if not Path(r["folder"]).is_dir():
            print(f"  skipped {r['name']}: folder {r['folder']} no longer exists")
            continue
        ingest(r["folder"], r["name"], r["commit_sha"], r["note"], r["targets_file"], r["started"])
    con = connect()
    with con:
        for n in notes:
            con.execute("INSERT INTO notes (created, star, run, candidate, tag, text) VALUES (?,?,?,?,?,?)",
                        (n["created"], n["star"], n["run"], n["candidate"], n["tag"], n["text"]))
        for t, rows in kept.items():
            for r in rows:
                con.execute(f"INSERT OR REPLACE INTO {t} ({', '.join(r)}) VALUES ({', '.join('?' * len(r))})",
                            tuple(r.values()))
    con.close()
    keep = BACKUPS / f"tess-before-rebuild-{datetime.datetime.now():%Y%m%d-%H%M%S}.db"
    keep.parent.mkdir(parents=True, exist_ok=True)
    os.replace(final, keep)
    os.replace(DB, final)
    DB = final
    print(f"rebuilt {final}: {len(runs)} runs, {len(notes)} notes (old copy: {keep})")
    print("restart `tessdb.py serve` if it is running")


def remove_run(name):
    """Delete a run and everything loaded from it (searches, candidates, flags,
    checks, injections). Stars no other run searched are removed too; notes stay."""
    con = connect()
    with con:
        n = con.execute("DELETE FROM runs WHERE name = ?", (name,)).rowcount
        if n:
            con.execute("DELETE FROM stars WHERE star NOT IN (SELECT star FROM searches) "
                        "AND star NOT IN (SELECT star FROM notes)")
    print(f"removed run {name}" if n else f"no run named {name}")


def add_note(star, text, tag=None, run=None, candidate=None):
    con = connect()
    with con:
        con.execute("INSERT INTO notes (created, star, run, candidate, tag, text) VALUES (?,?,?,?,?,?)",
                    (datetime.datetime.now().isoformat(timespec="seconds"), star, run, candidate, tag, text))
    print(f"note added for {star}")


def list_runs():
    con = connect()
    rows = con.execute("SELECT name, kind, started, substr(commit_sha, 1, 7), n_targets,"
                       " (SELECT count(*) FROM candidates c WHERE c.run_id = r.run_id),"
                       " (SELECT count(*) FROM candidates c WHERE c.run_id = r.run_id AND status = 'clean')"
                       " FROM runs r ORDER BY started").fetchall()
    print(f"{'run':28s} {'kind':17s} {'started':19s} {'commit':7s} {'stars':>5s} {'signals':>7s} {'clean':>5s}")
    for r in rows:
        print(f"{r[0]:28s} {r[1]:17s} {(r[2] or ''):19s} {(r[3] or '-'):7s} {r[4] or 0:5d} {r[5]:7d} {r[6]:5d}")


def backup():
    out = BACKUPS / f"tess-{datetime.date.today()}.db"
    out.parent.mkdir(parents=True, exist_ok=True)
    src = sqlite3.connect(DB)
    with sqlite3.connect(out) as dst:
        src.backup(dst)                          # consistent copy, safe while in use
    print(f"backed up to {out}")
    export_records()


def export_records():
    """Write RECORD_TABLES to RECORDS/<table>.csv (notes without their row ids, which
    a rebuild renumbers)."""
    RECORDS.mkdir(parents=True, exist_ok=True)
    con = connect()
    counts = []
    for t, order in RECORD_TABLES.items():
        cur = con.execute(f"SELECT * FROM {t} ORDER BY {order}")
        cols = [d[0] for d in cur.description]
        keep = [i for i, c in enumerate(cols) if c != "note_id"]
        rows = cur.fetchall()
        with open(RECORDS / f"{t}.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow([cols[i] for i in keep])
            w.writerows([[r[i] for i in keep] for r in rows])
        counts.append(f"{len(rows)} {t}")
    print(f"records written to {RECORDS}/: " + ", ".join(counts) + " (commit them to keep them in git)")


def restore_records(folder=None):
    """Load RECORDS/*.csv back into the database, e.g. into a new database after the
    old one was lost. Notes already present (same time, star, tag and text) are
    skipped; other rows replace any with the same key."""
    folder = Path(folder) if folder else RECORDS
    con = connect()
    with con:
        for t in RECORD_TABLES:
            f = folder / f"{t}.csv"
            if not f.exists():
                continue
            rows = read_csv(f)
            n = 0
            for r in rows:
                r = {k: (v if v != "" else None) for k, v in r.items()}
                if t == "notes":
                    if con.execute("SELECT 1 FROM notes WHERE created IS ? AND star IS ? AND tag IS ? AND text IS ?",
                                   (r["created"], r["star"], r["tag"], r["text"])).fetchone():
                        continue
                con.execute(f"INSERT OR REPLACE INTO {t} ({', '.join(r)}) VALUES ({', '.join('?' * len(r))})",
                            tuple(r.values()))
                n += 1
            print(f"  {t}: {n} of {len(rows)} row(s) loaded")


def serve(port):
    """Datasette on 127.0.0.1, with results/ mounted for plots (so runs that start
    after the server, e.g. a batch with rolling updates, have working plot links),
    plus each run folder elsewhere (scratch runs) under its own name."""
    con = connect()
    mounts = ["--static", f"results:{HERE / 'results'}"]
    for name, folder in con.execute("SELECT name, folder FROM runs"):
        if folder and Path(folder).is_dir():
            mounts += ["--static", f"{name}:{folder}"]
    meta = {"title": "PlanetFinder", "databases": {DB.stem: {
        "description": "Start with candidates_view: every signal with its run, worker, flags, "
                       "notes and a clickable plot_url. Plots of runs ingested after this "
                       "server started need a restart of `tessdb.py serve`.",
        "tables": {
        "candidates_view": {"description": "Every signal: flags, notes, a plot link (click plot_url), the transit fit's density check and the TRICERATOPS false-positive probability",
                            "facets": ["status", "run"]},
        "stars_searched": {"description": "Each star searched, in how many runs, and its signals"},
        "watchlist": {"description": "Candidates re-tested as new TESS sectors arrive (watchlist.py check): "
                      "verdict from the events in sectors released after the candidate was found"},
        "vsx_shortlist_view": {"description": "Variable-star signals for the weekly VSX session "
                               "(vsx_shortlist.py): group 'new' = not in VSX, 'published' = in the TESS Ten Thousand Catalog (2025ApJS..279...50K) at this period, 'period' = in VSX "
                               "without the TESS period. notes shows stars already queued or submitted",
                               "facets": ["group", "run"]},
        "catalog": {"description": "Whole-sky target catalogue: score = expected SNR of a "
                               "2 R_earth, 10 d planet"},
        "runs": {"label_column": "name"},
        "candidates": {"facets": ["status", "run_id"]},
        "candidate_flags": {"facets": ["flag", "source"]},
        "notes": {"facets": ["tag"]}}}}}
    meta_path = DATASETTE / "metadata.json"
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    # wide text columns (flags, notes) otherwise push the table far off-screen, with its
    # scroll bar out of sight below the rows: wrap cell text and cap column widths
    assets = DATASETTE / "assets"
    assets.mkdir(exist_ok=True)
    (assets / "tess.css").write_text(
        ".rows-and-columns td { white-space: normal; word-break: break-word; max-width: 22em; "
        "vertical-align: top; }\n"
        ".rows-and-columns td.col-plot_url { max-width: 30em; font-size: 0.85em; }\n"
        ".table-wrapper { overflow-x: auto; }\n")
    meta["extra_css_urls"] = ["/tess-assets/tess.css"]
    mounts += ["--static", f"tess-assets:{assets}"]
    meta_path.write_text(json.dumps(meta))
    exe = Path(sys.executable).with_name("datasette")
    cmd = [str(exe), "serve", str(DB), "--host", "127.0.0.1", "--port", str(port),
           "--metadata", str(meta_path), "--setting", "sql_time_limit_ms", "5000",
           "--setting", "truncate_cells_html", "80", *mounts]
    print(f"Datasette on http://localhost:{port}  (Ctrl-C to stop)")
    subprocess.run(cmd)


def clean_cache(dry_run=False):
    """Delete cached light curves (lightkurve's mastDownload folder) of stars whose
    every search ended with no signal, and that have no notes. find_planets.py does
    this itself for new searches; this tidies up stars searched before it did.
    Stars with any signal, in any run, or with notes keep their data."""
    from lightkurve.config import get_cache_dir
    con = connect()
    status = {}
    for star, st in con.execute("SELECT star, status FROM searches"):
        status.setdefault(star, set()).add(st)
    keep = ({r[0] for r in con.execute("SELECT DISTINCT star FROM candidates")}
            | {r[0] for r in con.execute("SELECT DISTINCT star FROM notes")})
    drop = set()
    for star, sts in status.items():
        m = re.fullmatch(r"TIC (\d+)", star)
        if m and star not in keep and sts <= {"ok_0_candidates", "no_good_data"}:
            drop.add(int(m.group(1)))
    root = Path(get_cache_dir()).resolve() / "mastDownload"
    files = [f for f in root.rglob("*.fits")
             if any(int(x) in drop for x in re.findall(r"(?<!\d)\d{16}(?!\d)", f.name))]
    size = sum(f.stat().st_size for f in files)
    stars = {int(x) for f in files for x in re.findall(r"(?<!\d)\d{16}(?!\d)", f.name)} & drop
    print(f"{len(drop)} stars searched with no signal and no notes; {len(stars)} have cached "
          f"light curves: {len(files)} files, {size / 1e9:.2f} GB")
    if dry_run:
        print("dry run: nothing deleted")
        return
    for f in files:
        f.unlink()
    for d in sorted({f.parent for f in files}, key=lambda x: -len(x.parts)):
        while d != root and d.is_dir() and not any(d.iterdir()):
            d.rmdir()
            d = d.parent
    print(f"deleted {len(files)} files ({size / 1e9:.2f} GB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("ingest", help="load a finished results folder")
    a.add_argument("folder")
    a.add_argument("--name", help="run name (default: folder name)")
    a.add_argument("--commit", help="code version, for folders without run_info.json")
    a.add_argument("--note", help="free-text note about the run")
    a.add_argument("--targets", help="target list, for folders without run_info.json")
    a.add_argument("--started", help="start time (YYYY-MM-DD HH:MM:SS), if not in the folder")
    n = sub.add_parser("note", help="record a note against a star")
    n.add_argument("star")
    n.add_argument("text")
    n.add_argument("--tag")
    n.add_argument("--run")
    n.add_argument("--candidate", type=int)
    sub.add_parser("list", help="list runs")
    rm = sub.add_parser("remove", help="delete a run and everything loaded from it")
    rm.add_argument("name")
    sub.add_parser("backup", help="dated copy in db/backups/, and records/*.csv")
    rr = sub.add_parser("restore-records", help="load records/*.csv back into the database")
    rr.add_argument("folder", nargs="?", help="default: records/")
    sub.add_parser("rebuild", help="recreate with the current schema (re-ingests every run)")
    cc = sub.add_parser("clean-cache", help="delete cached light curves of stars with no signal")
    cc.add_argument("--dry-run", action="store_true", help="only report what would be deleted")
    s = sub.add_parser("serve", help="browse with Datasette")
    s.add_argument("--port", type=int, default=PORT)
    args = ap.parse_args()
    if args.cmd == "ingest":
        ingest(args.folder, args.name, args.commit, args.note, args.targets, args.started)
    elif args.cmd == "note":
        add_note(args.star, args.text, args.tag, args.run, args.candidate)
    elif args.cmd == "list":
        list_runs()
    elif args.cmd == "remove":
        remove_run(args.name)
    elif args.cmd == "backup":
        backup()
    elif args.cmd == "restore-records":
        restore_records(args.folder)
    elif args.cmd == "rebuild":
        rebuild()
    elif args.cmd == "clean-cache":
        clean_cache(args.dry_run)
    elif args.cmd == "serve":
        serve(args.port)


if __name__ == "__main__":
    main()
