#!/usr/bin/env python3
"""
prefetch.py - fetch a batch's MAST data ahead of its search.

Usage:
    python prefetch.py targets/<name>.txt results/<name> [--workers 4]

For each target, does exactly what find_planets.py does before it searches: the
MAST product search, the choice of light curves (choose_products), the downloads
into lightkurve's cache, and the TIC lookup of radius and mass. Each star's result
goes into <out>/prefetch.jsonl (one JSON line per star: the files, the sources
label, the stellar parameters and when it was fetched). A search run on <out>
then reads a prefetched star entirely from the cache, without contacting MAST
(run_parallel.py points the workers at the manifest), so a MAST outage while the
batch runs no longer stalls it. Stars missing from the manifest, or whose files
have since gone, are searched the usual way.

Run it whenever MAST is healthy, for example when the batch is created. Rerunning
fetches only what's missing or failed, so it can be repeated until complete. The
light curves take about 6 MB a star until the search deletes those of stars with
no signal. A manifest more than STALE_DAYS old may miss newly released sectors;
rerun with --refresh to fetch again. While it runs it holds <out>/.prefetch.lock: a
second prefetch on the same folder is refused (exit 4), and a search started meanwhile
waits for it to finish.
"""
import argparse
import fcntl
import json
import logging
import sys
import time
import warnings
from datetime import datetime, timezone
from multiprocessing import Pool
from pathlib import Path

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
STALE_DAYS = 14


def read_manifest(path):
    """{star: record}, the last record per star."""
    out = {}
    if Path(path).exists():
        for line in Path(path).read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                out[r["star"]] = r
    return out


def usable(rec):
    """A record the search can use: fetched, with all its files still present."""
    return rec.get("status") == "ok" and all(Path(f).is_file() for f in rec.get("files", []))


def fetch_one(name):
    import find_planets as fp
    import lightkurve as lk
    from mastutil import retry
    logging.disable(logging.WARNING)
    t0 = time.time()
    rec = dict(star=name, fetched=datetime.now(timezone.utc).isoformat(timespec="seconds"))
    try:
        sr = retry(lambda: lk.search_lightcurve(name, author=fp.AUTHORS), f"{name}: MAST search")
        if len(sr) == 0:
            return dict(rec, status="no_data", files=[], seconds=round(time.time() - t0))
        rows, sources = fp.choose_products(sr)
        if not rows:
            return dict(rec, status="no_usable_data", files=[], seconds=round(time.time() - t0))
        lcc = fp.download(sr[rows])
        if lcc is None or len(lcc) == 0:
            return dict(rec, status="error: download_failed", seconds=round(time.time() - t0))
        files = [str(Path(lc.meta["FILENAME"]).resolve()) for lc in lcc]
        star = retry(lambda: fp.star_params(name, raise_errors=True), f"{name}: TIC lookup")
        return dict(rec, status="ok", files=files, sources=sources, star_params=star,
                    seconds=round(time.time() - t0))
    except Exception as e:
        return dict(rec, status=f"error: {type(e).__name__}: {e}"[:200], seconds=round(time.time() - t0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("targets")
    ap.add_argument("out", help="the batch's results folder (results/<name>)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--refresh", action="store_true", help="fetch every star again")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    import find_planets as fp
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    lock = open(out / ".prefetch.lock", "w")       # held until exit; run_parallel.py waits on it
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        logging.error(f"another prefetch is already running on {out}")
        return 4
    path = out / "prefetch.jsonl"
    have = {} if args.refresh else read_manifest(path)
    targets = fp.read_targets(args.targets)
    todo = [t for t in targets if not (t in have and (usable(have[t]) or have[t]["status"] in
                                                      ("no_data", "no_usable_data")))]
    logging.info(f"{len(targets)} targets: {len(targets) - len(todo)} already prefetched, "
                 f"{len(todo)} to fetch with {args.workers} workers -> {path}")
    n_ok = n_err = 0
    with open(path, "a") as fh, Pool(args.workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(fetch_one, todo), 1):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if rec["status"].startswith("error"):
                n_err += 1
                logging.warning(f"[{i}/{len(todo)}] {rec['star']}: {rec['status']}")
            else:
                n_ok += 1
            if i % 25 == 0 or i == len(todo):
                logging.info(f"[{i}/{len(todo)}] fetched {n_ok}, failed {n_err}")
    final = read_manifest(path)
    ready = sum(1 for t in targets if t in final and (usable(final[t]) or
                                                      final[t]["status"] in ("no_data", "no_usable_data")))
    logging.info(f"prefetch: {ready}/{len(targets)} targets ready"
                 + (f"; {len(targets) - ready} still to fetch (rerun to retry)" if ready < len(targets) else ""))
    return 0 if ready == len(targets) else 3


if __name__ == "__main__":
    sys.exit(main())
