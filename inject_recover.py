#!/usr/bin/env python3
"""
inject_recover.py - injection-recovery test of find_planets.py's search.

Usage:
    python inject_recover.py targets/ffi.txt [--n 500] [--workers 8] [--out results/injection]
    python inject_recover.py --summary results/injection
    python inject_recover.py targets/aries.txt --invert --out results/falsealarm_aries
    python inject_recover.py targets/controls.txt --rp 2.3 2.7 --period 10 25 --n 60 --out ...

Each job injects one synthetic transit into a target star's real (cleaned) TESS
light curves and runs the unchanged find_planets search on the result:
  - planet radius log-uniform in RP_RANGE (Earth radii), period log-uniform in
    P_RANGE (days), impact parameter uniform in [0, B_MAX], random epoch;
  - the transit is a trapezoid whose depth (Rp/R*)^2 and durations follow from
    the star's TIC radius and mass for a circular orbit;
  - recovered = a logged signal within REC_P_TOL of the injected period whose
    epoch lands within one duration of an injected transit; a signal at 2P or
    P/2 is counted separately as an alias.
The expected SNR (depth over the per-point scatter, summed over in-transit
points) is recorded, so completeness can be read against SNR as well as against
radius and period. Jobs come from a fixed seed and results are appended as they
finish: rerun the same command to resume. --summary prints completeness tables.
--rp and --period override RP_RANGE and P_RANGE.

--invert is a false-alarm test instead: each star's light curves are flipped
(2 - flux), so real transits become brightenings the search ignores while noise
and spot variability, roughly symmetric, still make dips. Every signal found is
then a false alarm; their statuses go to falsealarms.csv.
"""
import argparse
import csv
import io
import logging
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")       # one thread per worker process
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np                                    # noqa: E402

RP_RANGE = (1.0, 8.0)       # Earth radii (faint field stars: smaller is mostly undetectable)
P_RANGE = (0.5, 30.0)       # days
B_MAX = 0.9
MIN_TRANSITS = 3            # injections with fewer transits in the data are redrawn
REC_P_TOL = 0.005           # fractional period match for a recovery
SEED = 20261011
RESULT_FIELDS = ["job", "star", "rstar", "mstar", "p_inj", "t0_inj", "rp_inj", "b_inj",
                 "depth_inj_ppm", "t14_h", "n_tr_inj", "snr_exp", "recovered", "alias",
                 "rec_status", "rec_flags", "rec_p", "rec_depth_ppm", "rec_sde", "rec_sde_margin",
                 "rec_model_shift", "n_signals", "minutes"]

_cache = {}                 # per worker process: star -> (lcs, star_params, sectors, sources)


def read_targets(path):
    out = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def load_star(name):
    """Cleaned light curves etc. for a star, cached in this worker process."""
    import find_planets as fp
    import lightkurve as lk
    if name not in _cache:
        sr = lk.search_lightcurve(name, author=fp.AUTHORS)
        rows, sources = fp.choose_products(sr) if len(sr) else ([], "")
        lcc = fp.download(sr[rows]) if rows else None
        lcs = [c for c in (fp.clean(lc) for lc in (lcc or [])) if c is not None]
        star = fp.star_params(name)
        sectors = [(f"{c.meta.get('SECTOR')}{'(QLP)' if c.meta.get('AUTHOR') == 'QLP' else ''}",
                    float(fp.arr(c.time).min()), float(fp.arr(c.time).max())) for c in lcs]
        _cache[name] = (lcs, star, sectors, sources)
    return _cache[name]


def trapezoid(t, P, T0, depth, t14, t23):
    """Fractional flux drop of a trapezoidal transit at times t."""
    dt = np.abs((t - T0 + 0.5 * P) % P - 0.5 * P)
    half14, half23 = t14 / 2, t23 / 2
    drop = np.zeros_like(t)
    drop[dt <= half23] = depth
    ramp = (dt > half23) & (dt < half14)
    drop[ramp] = depth * (half14 - dt[ramp]) / max(half14 - half23, 1e-9)
    return drop


def make_jobs(stars, n, seed=SEED, rp_range=RP_RANGE, p_range=P_RANGE):
    """Deterministic list of (job, star, uniform draws); physical parameters are
    worked out per star in the worker, where its radius and data are known."""
    rng = np.random.default_rng(seed)
    return [dict(job=i, star=stars[rng.integers(len(stars))], rp_range=tuple(rp_range),
                 p_range=tuple(p_range),
                 u_rp=rng.random(), u_p=rng.random(), u_b=rng.random(), u_t0=rng.random())
            for i in range(n)]


FA_FIELDS = ["job", "star", "n_signals", "n_clean", "n_flagged", "signals", "minutes"]


def run_invert(job):
    """False-alarm test: search the star's flipped light curves; report every signal."""
    import lightkurve as lk
    from astropy.time import Time
    import find_planets as fp
    logging.disable(logging.CRITICAL)
    fp.save_plot = lambda *a, **k: None
    t_start = time.time()
    name = job["star"]
    lcs, star, sectors, sources = load_star(name)
    if not lcs:
        return dict(job=job["job"], star=name, signals="no data")
    flipped = [lk.LightCurve(time=Time(fp.arr(c.time), format=c.time.format, scale=c.time.scale),
                             flux=2 - fp.arr(c.flux), flux_err=fp.arr(c.flux_err), meta=dict(c.meta))
               for c in lcs]
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fp.CSV_FIELDS)
    outdir = Path(os.environ.get("TMPDIR", "/tmp")) / f"invert_{os.getpid()}"
    fp.search_star(flipped, name, len(flipped), outdir, writer, buf, sources, star, sectors)
    rows = list(csv.DictReader(io.StringIO(buf.getvalue()), fieldnames=fp.CSV_FIELDS))
    rows = [r for r in rows if r["status"] != "same_period"]
    return dict(job=job["job"], star=name, n_signals=len(rows),
                n_clean=sum(r["status"] == "clean" for r in rows),
                n_flagged=sum(r["status"] == "flagged" for r in rows),
                signals="; ".join(f"P={float(r['period_d']):.4f} d {float(r['depth_ppm']):.0f} ppm "
                                  f"SNR {float(r['snr']):.1f} {r['status']} {r['flags']}" for r in rows),
                minutes=round((time.time() - t_start) / 60, 2))


def run_job(job):
    if job.get("invert"):
        return run_invert(job)
    import lightkurve as lk
    from astropy.time import Time
    import find_planets as fp
    logging.disable(logging.CRITICAL)                  # the search logs per pass; keep quiet
    fp.save_plot = lambda *a, **k: None                # plots are not needed here
    t_start = time.time()
    name = job["star"]
    lcs, star, sectors, sources = load_star(name)
    base = dict(job=job["job"], star=name)
    if not lcs or not star:
        return dict(base, recovered="skip", rec_flags="no data or no TIC radius/mass")
    t_all = np.concatenate([fp.arr(c.time) for c in lcs])
    span_lo, span_hi = t_all.min(), t_all.max()
    rng = np.random.default_rng(SEED + 7919 * job["job"])
    # physical parameters; redraw period/epoch until at least MIN_TRANSITS fall on data
    RP, PR = job.get("rp_range", RP_RANGE), job.get("p_range", P_RANGE)
    rp = RP[0] * (RP[1] / RP[0]) ** job["u_rp"]
    b = B_MAX * job["u_b"]
    k = rp / 109.1 / star["rad"]
    u_p, u_t0 = job["u_p"], job["u_t0"]
    for _ in range(50):
        P = PR[0] * (PR[1] / PR[0]) ** u_p
        a = 215.032 * (star["mass"] * (P / 365.25) ** 2) ** (1 / 3) / star["rad"]   # a/R*
        x14 = (1 + k) ** 2 - b ** 2
        T0 = span_lo + u_t0 * P
        if x14 > 0 and a > 1 + k:
            t14 = P / np.pi * np.arcsin(min(1.0, np.sqrt(x14) / a))
            x23 = (1 - k) ** 2 - b ** 2
            t23 = P / np.pi * np.arcsin(min(1.0, np.sqrt(x23) / a)) if x23 > 0 else 0.0
            n_tr = len({int(np.round((tt - T0) / P)) for tt in
                        t_all[np.abs((t_all - T0 + 0.5 * P) % P - 0.5 * P) < 0.5 * t14]})
            if n_tr >= MIN_TRANSITS:
                break
        u_p, u_t0 = rng.random(), rng.random()
    else:
        return dict(base, recovered="skip", rec_flags="no geometry with enough transits")
    depth = k ** 2
    # inject; expected SNR from each sector's point-to-point scatter
    inj, snr2 = [], 0.0
    for c in lcs:
        t, f, e = fp.arr(c.time), fp.arr(c.flux), fp.arr(c.flux_err)
        drop = trapezoid(t, P, T0, depth, t14, t23)
        sig = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
        snr2 += np.sum(drop ** 2) / sig ** 2
        inj.append(lk.LightCurve(time=Time(t, format=c.time.format, scale=c.time.scale),
                                 flux=f * (1 - drop), flux_err=e, meta=dict(c.meta)))
    # search exactly as find_planets does
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fp.CSV_FIELDS)
    outdir = Path(os.environ.get("TMPDIR", "/tmp")) / f"inject_{os.getpid()}"
    fp.search_star(inj, name, len(inj), outdir, writer, buf, sources, star, sectors)
    rows = list(csv.DictReader(io.StringIO(buf.getvalue()), fieldnames=fp.CSV_FIELDS))
    recovered, alias, best = False, False, None
    for r in rows:
        p, t0r = float(r["period_d"]), float(r["t0_btjd"])
        d = max(float(r["duration_h"]) / 24, t14)
        phase_ok = abs((t0r - T0 + 0.5 * P) % P - 0.5 * P) < d
        if abs(p / P - 1) < REC_P_TOL and phase_ok:
            recovered, best = True, r
            break
        if any(abs(p / (h * P) - 1) < REC_P_TOL for h in (2.0, 0.5)) and not alias:
            alias, best = True, best or r
    return dict(base, rstar=round(star["rad"], 3), mstar=round(star["mass"], 3),
                p_inj=round(P, 5), t0_inj=round(T0, 5), rp_inj=round(rp, 3), b_inj=round(b, 3),
                depth_inj_ppm=round(depth * 1e6, 1), t14_h=round(t14 * 24, 2), n_tr_inj=n_tr,
                snr_exp=round(float(np.sqrt(snr2)), 1), recovered=recovered, alias=alias,
                rec_status=best["status"] if best else "", rec_flags=best["flags"] if best else "",
                rec_p=best["period_d"] if best else "",
                rec_depth_ppm=best["depth_ppm"] if best else "", rec_sde=best["sde"] if best else "",
                rec_sde_margin=best["sde_margin"] if best else "",
                rec_model_shift=best["model_shift"] if best else "", n_signals=len(rows),
                minutes=round((time.time() - t_start) / 60, 2))


def summary(outdir):
    import pandas as pd
    d = pd.read_csv(Path(outdir) / "injections.csv")
    d = d[d.recovered.isin(["True", "False", True, False])].copy()
    d["rec"] = d.recovered.astype(str) == "True"
    d["clean"] = d.rec & (d.rec_status == "clean")
    print(f"{len(d)} injections; recovered {d.rec.mean() * 100:.0f}% "
          f"(as clean {d.clean.mean() * 100:.0f}%), aliases {(d.alias.astype(str) == 'True').mean() * 100:.0f}%")
    d["snr_bin"] = pd.cut(d.snr_exp, [0, 5, 7, 10, 15, 20, 30, 50, 1e9])
    print("\nby expected SNR:")
    print(d.groupby("snr_bin", observed=True).agg(n=("rec", "size"), recovered=("rec", "mean"),
                                                  clean=("clean", "mean")).round(2).to_string())
    d["rp_bin"] = pd.cut(d.rp_inj, [1, 1.5, 2, 3, 4, 6, 8])
    d["p_bin"] = pd.cut(d.p_inj, [0.5, 1, 3, 10, 30])
    print("\nrecovered fraction by radius (rows, Earth radii) and period (columns, days):")
    print(d.pivot_table(index="rp_bin", columns="p_bin", values="rec", aggfunc="mean",
                        observed=True).round(2).to_string())
    print("\ncounts:")
    print(d.pivot_table(index="rp_bin", columns="p_bin", values="rec", aggfunc="size",
                        observed=True).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("targets", nargs="?", help="target list (stars to inject into)")
    ap.add_argument("--n", type=int, default=500, help="number of injections (default 500)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default="results/injection")
    ap.add_argument("--rp", nargs=2, type=float, metavar=("MIN", "MAX"), help="radius range, Earth radii")
    ap.add_argument("--period", nargs=2, type=float, metavar=("MIN", "MAX"), help="period range, days")
    ap.add_argument("--invert", action="store_true", help="false-alarm test on flipped light curves")
    ap.add_argument("--summary", metavar="DIR", help="print completeness tables and exit")
    args = ap.parse_args()
    if args.summary:
        summary(args.summary)
        return
    if not args.targets:
        sys.exit("give a target list, or --summary DIR")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / ("falsealarms.csv" if args.invert else "injections.csv")
    fields = FA_FIELDS if args.invert else RESULT_FIELDS
    done = set()
    if path.exists():
        with open(path) as f:
            rd = csv.DictReader(f)
            done = {int(r["job"]) for r in rd}
            fields = [c for c in rd.fieldnames]       # resume with the file's own columns
    if args.invert:                                   # one job per star
        jobs = [dict(job=i, star=st, invert=True) for i, st in enumerate(read_targets(args.targets))
                if i not in done]
    else:
        jobs = [j for j in make_jobs(read_targets(args.targets), args.n,
                                     rp_range=args.rp or RP_RANGE, p_range=args.period or P_RANGE)
                if j["job"] not in done]
    # group jobs by star so each worker reuses its cached light curves
    jobs.sort(key=lambda j: (j["star"], j["job"]))
    print(f"{len(done)} done, {len(jobs)} to run with {args.workers} workers -> {path}", flush=True)
    import multiprocessing as mp
    new = not path.exists()
    with open(path, "a", newline="") as fh, mp.get_context("fork").Pool(args.workers) as pool:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if new:
            w.writeheader()
        for i, res in enumerate(pool.imap_unordered(run_job, jobs, chunksize=2), 1):
            w.writerow({k: res.get(k, "") for k in fields})
            fh.flush()
            if i % 10 == 0 or i == len(jobs):
                print(f"{i}/{len(jobs)} {time.strftime('%H:%M')}", flush=True)


if __name__ == "__main__":
    main()
