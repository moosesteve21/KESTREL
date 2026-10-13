#!/usr/bin/env python3
"""
watchlist.py - re-test candidates as new TESS sectors are released.

Usage:
    python watchlist.py add "TIC 123456789" --run b011 --candidate 1 --reason "doubtful: S28/S68 missing"
    python watchlist.py add "TIC 123" --period 7.1 --t0 1331.7 --duration-h 2.3 --depth-ppm 300 \\
                            --sectors 1,28 --reason "..."
    python watchlist.py list
    python watchlist.py check [STAR ...]          # all watched stars by default
    python watchlist.py remove "TIC 123"

A candidate found in some sectors makes a prediction for every later one: a dip of
known depth and duration at a known time. `check` asks MAST which sectors exist now,
downloads the ones not yet measured (same products as the search: SPOC, then
TESS-SPOC, then QLP) and measures the depth of each transit at its predicted time,
with the shape held fixed and a local quadratic baseline. Each sector's noise is
also measured at the same duration at times away from the transits, and the larger
of that and the formal error is used, so red noise doesn't inflate the result.

Sectors the candidate was found in are "baseline": measured the same way (to show
the method sees the signal there) but left out of the verdict. The verdict uses only
later sectors, an independent test:
    supports      new depth > 0 at >= 3 sigma and within 3 sigma of the prediction
    contradicts   new depth below the prediction by >= 3 sigma
    inconclusive  neither yet: the report says roughly how many more transits would settle it
    no new data   no transit in any later sector yet
The report also lists the sectors tess-point says will observe the star next.

Results go into tess.db (watchlist, watch_sectors, watch_events) and a dated report
with event plots in db/watchlist/. Run `check` after new sectors come out (TESS
light curves reach MAST a month or two after each ~27-day sector; @MAST_News
announces ingests).
"""
import argparse
import datetime
import json
import logging
import re
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tessdb  # noqa: E402

REPORTS = tessdb.DB.parent / "watchlist"
TAU_FRAC = 0.15             # ingress/egress as a fraction of T14 (the shape is held fixed)
WINDOW_T14 = 3.0            # each event's baseline window: +- this many T14 (at least 0.75 d)
MIN_IN, MIN_OUT = 5, 15     # points needed in transit, and on each side out of transit
NULL_STEP_T14 = 1.5         # null-test spacing (in T14), kept this far from real transits
SIGMA = 3.0


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def slug(star):
    return re.sub(r"[^\w.-]+", "_", star)


def parse_sectors(sources):
    """'SPOC/120s: 1; TESS-SPOC/600s: 28; QLP/200s: 92,95' -> [1, 28, 92, 95]."""
    return sorted({int(x) for part in str(sources or "").split(";") if ":" in part
                   for x in part.split(":", 1)[1].split(",") if x.strip().isdigit()})


# ---------------------------------------------------------------------- add --
def add(con, star, args):
    rec = dict(star=star, reason=args.reason or "", added=now())
    if args.run:
        row = con.execute("SELECT c.period_d, c.t0_btjd, c.duration_h, c.depth_ppm, c.data, r.run_id "
                          "FROM candidates c JOIN runs r USING (run_id) WHERE r.name = ? AND c.star = ? "
                          "AND c.candidate = ?", (args.run, star, args.candidate)).fetchone()
        if not row:
            sys.exit(f"no candidate {args.candidate} of {star} in run {args.run}")
        P, T0, dur, depth, data, run_id = row
        data = json.loads(data or "{}")
        rec.update(period=P, t0=T0, duration_h=dur, depth_ppm=depth, period_err=None, t0_err=None,
                   baseline_sectors=",".join(map(str, parse_sectors(data.get("sources")))),
                   source=f"{args.run} candidate {args.candidate} (search)")
        fit = con.execute("SELECT data FROM transit_fits WHERE run_id = ? AND star = ? AND candidate = ?",
                          (run_id, star, args.candidate)).fetchone()
        if fit:                                         # the MCMC fit's values are better
            f = {k: float(v) for k, v in json.loads(fit[0]).items()
                 if re.fullmatch(r"(period|t0|t14_h|depth_ppm)(_lo|_hi)?", k) and v not in ("", None)}
            if "period" in f:
                rec.update(period=f["period"], t0=f["t0"], duration_h=f.get("t14_h", dur),
                           depth_ppm=f.get("depth_ppm", depth),
                           period_err=(f["period_hi"] - f["period_lo"]) / 2 if "period_hi" in f else None,
                           t0_err=(f["t0_hi"] - f["t0_lo"]) / 2 if "t0_hi" in f else None,
                           source=f"{args.run} candidate {args.candidate} (transit fit)")
    else:
        if None in (args.period, args.t0, args.duration_h, args.depth_ppm):
            sys.exit("give --run/--candidate, or all of --period --t0 --duration-h --depth-ppm")
        rec.update(period=args.period, t0=args.t0, duration_h=args.duration_h, depth_ppm=args.depth_ppm,
                   period_err=args.period_err, t0_err=args.t0_err, source="given")
    if args.sectors is not None:
        rec["baseline_sectors"] = args.sectors
    rec.setdefault("baseline_sectors", "")
    with con:
        con.execute("DELETE FROM watchlist WHERE star = ?", (star,))
        con.execute(f"INSERT INTO watchlist ({', '.join(rec)}) VALUES ({', '.join('?' * len(rec))})",
                    tuple(rec.values()))
    print(f"watching {star}: P = {rec['period']:.6f} d, T0 = {rec['t0']:.4f}, {rec['duration_h']:.2f} h, "
          f"{rec['depth_ppm']:.0f} ppm; baseline sectors {rec['baseline_sectors'] or 'none'} ({rec['source']})")


# -------------------------------------------------------------- measurement --
def shape(dt, t14, tau):
    a, half = np.abs(dt), t14 / 2
    return np.clip((half - a) / tau, 0, 1)


def measure(t, y, e, tc, t14, cad):
    """Depth (fraction) and formal error of a fixed-shape dip at tc, with a local
    quadratic baseline; None if the event isn't covered."""
    win = max(0.75, WINDOW_T14 * t14)
    sel = np.abs(t - tc) < win
    if not sel.any():
        return None
    dt = t[sel] - tc
    if (np.sum(np.abs(dt) < t14 / 2) < MIN_IN or np.sum(dt < -0.75 * t14) < MIN_OUT
            or np.sum(dt > 0.75 * t14) < MIN_OUT):
        return None
    sub = (np.arange(5) - 2) / 5 * cad                     # integrate over the exposure
    s = np.mean([shape(dt + o, t14, TAU_FRAC * t14) for o in sub], axis=0)
    X = np.column_stack([np.ones_like(dt), dt, dt ** 2, -s])
    w = 1 / e[sel]
    c, *_ = np.linalg.lstsq(X * w[:, None], y[sel] * w, rcond=None)
    r = (y[sel] - X @ c) * w
    cov = np.linalg.inv((X * w[:, None]).T @ (X * w[:, None])) * max(1.0, np.sum(r ** 2) / max(1, len(r) - 4))
    return c[3] / c[0], np.sqrt(cov[3, 3]) / c[0]


def sector_events(lc, w, sector, product):
    import find_planets as fp
    t, y, e = fp.arr(lc.time), fp.arr(lc.flux), fp.arr(lc.flux_err)
    cad = np.median(np.diff(t))
    t14 = w["duration_h"] / 24
    P, T0 = w["period"], w["t0"]
    epochs = np.arange(np.ceil((t.min() - T0) / P - 0.1), np.floor((t.max() - T0) / P + 0.1) + 1)
    # null test: same measurement at times away from any transit, for this sector's real noise
    real = T0 + epochs * P
    fake = np.arange(t.min() + 1, t.max() - 1, NULL_STEP_T14 * t14)
    fake = fake[np.min(np.abs(fake[:, None] - real[None, :]), axis=1) > 2 * t14] if len(real) else fake
    null = [m[0] for m in (measure(t, y, e, x, t14, cad) for x in fake) if m]
    sig_null = np.std(null) if len(null) >= 8 else 0.0
    out = []
    for k in epochs:
        tc = T0 + k * P
        m = measure(t, y, e, tc, t14, cad)
        if m is None:
            continue
        err = max(m[1], sig_null)
        out.append(dict(sector=sector, epoch=int(k), tc=tc, depth_ppm=m[0] * 1e6, err_ppm=err * 1e6,
                        snr=m[0] / err, product=product, t=t, y=y))
    return out, sig_null * 1e6, len(null)


def available(star):
    """{sector: (row index, product label)} of the products the search would use, and the result."""
    import find_planets as fp
    import lightkurve as lk
    from mastutil import retry, set_mast_timeout
    set_mast_timeout()
    sr = retry(lambda: lk.search_lightcurve(star, author=fp.AUTHORS), f"{star}: MAST search")
    rows, _ = fp.choose_products(sr)
    tab = sr.table
    return {int(tab["sequence_number"][i]): (i, f"{tab['author'][i]}/{int(tab['exptime'][i])}s")
            for i in rows}, sr


def upcoming(con, star, latest):
    """Sectors after `latest` that tess-point expects to observe the star."""
    try:
        from tess_stars2px import tess_stars2px_function_entry
        r = con.execute("SELECT ra, dec FROM stars WHERE star = ?", (star,)).fetchone()
        if not r or r[0] is None:
            return ""
        tic = int(re.sub(r"\D", "", star))
        _, _, _, sec, *_ = tess_stars2px_function_entry([tic], [r[0]], [r[1]])
        nxt = sorted({int(s) for s in sec if latest < s < 1000})
        return ",".join(map(str, nxt[:6]))
    except Exception as e:                                  # informative only
        logging.warning(f"{star}: tess-point failed ({e})")
        return ""


def wmean(ev):
    if not ev:
        return None
    d = np.array([x["depth_ppm"] for x in ev]); s = np.array([x["err_ppm"] for x in ev])
    w = 1 / s ** 2
    return float(np.sum(d * w) / w.sum()), float(1 / np.sqrt(w.sum())), float(np.median(s))


def verdict(new, expected):
    """(verdict, text) from the post-discovery events."""
    m = wmean(new)
    if m is None:
        return "no new data", "no transit covered in a later sector yet"
    mean, err, per_event = m
    txt = f"{len(new)} new transit(s): {mean:.0f} +- {err:.0f} ppm (predicted {expected:.0f})"
    if mean / err >= SIGMA and abs(mean - expected) < SIGMA * err:
        return "supports", txt
    if (expected - mean) / err >= SIGMA:
        return "contradicts", txt
    # transits needed (at this per-event noise) to tell `expected` from 0 at 3 sigma
    need = int(np.ceil((SIGMA * per_event / expected) ** 2)) if expected > 0 else 0
    return "inconclusive", txt + f"; about {max(need - len(new), 1)} more transit(s) needed at this noise"


# -------------------------------------------------------------------- check --
def check_star(con, w, report, plots):
    import lightkurve as lk
    import find_planets as fp
    star = w["star"]
    baseline = {int(s) for s in str(w["baseline_sectors"] or "").split(",") if s.strip()}
    done = {r[0] for r in con.execute("SELECT sector FROM watch_sectors WHERE star = ?", (star,))}
    have, sr = available(star)
    todo = sorted(set(have) - done)
    events = []
    for sec in todo:
        i, product = have[sec]
        try:
            lc = fp.clean(sr[i].download(quality_bitmask="default"))
        except Exception as e:
            logging.warning(f"{star} S{sec}: download failed ({type(e).__name__}); will retry next check")
            continue
        if lc is None:
            ev, sig_null, n_null = [], 0, 0
        else:
            ev, sig_null, n_null = sector_events(lc, w, sec, product)
        isbase = int(sec in baseline)
        with con:
            con.execute("INSERT OR REPLACE INTO watch_sectors VALUES (?,?,?,?,?,?)",
                        (star, sec, product, isbase, len(ev), now()))
            for x in ev:
                con.execute("INSERT OR REPLACE INTO watch_events VALUES (?,?,?,?,?,?,?,?,?,?)",
                            (star, sec, x["epoch"], x["tc"], x["depth_ppm"], x["err_ppm"], x["snr"],
                             isbase, product, now()))
        events += [dict(x, baseline=isbase) for x in ev]
        logging.info(f"  {star} S{sec} ({product}{', baseline' if isbase else ''}): {len(ev)} transit(s)"
                     + (f", null-test noise {sig_null:.0f} ppm from {n_null} positions" if n_null else ""))
    rows = [dict(zip(("sector", "epoch", "depth_ppm", "err_ppm", "baseline", "product"), r)) for r in con.execute(
        "SELECT sector, epoch, depth_ppm, err_ppm, baseline, product FROM watch_events WHERE star = ? "
        "ORDER BY sector, epoch", (star,))]
    new = [r for r in rows if not r["baseline"]]
    base = [r for r in rows if r["baseline"]]
    v, txt = verdict(new, w["depth_ppm"])
    latest = max(have) if have else 0
    nxt = upcoming(con, star, max([s for s in have if s < 1000], default=0))
    with con:
        con.execute("UPDATE watchlist SET last_checked = ?, verdict = ?, summary = ?, next_sectors = ? "
                    "WHERE star = ?", (now(), v, txt, nxt or None, star))
    # report
    bm = wmean(base)
    report.append(f"## {star}: {v}\n")
    report.append(f"- Candidate: P = {w['period']:.6f} d, T0 = {w['t0']:.4f} BTJD, {w['duration_h']:.2f} h, "
                  f"{w['depth_ppm']:.0f} ppm ({w['source']}). Watched because: {w['reason'] or '-'}")
    report.append(f"- New data since discovery: {txt}")
    if bm:
        report.append(f"- Baseline sectors, same method: {len(base)} transit(s), {bm[0]:.0f} +- {bm[1]:.0f} ppm")
    report.append(f"- Sectors measured this check: {', '.join(f'S{s}' for s in todo) or 'none new'}; "
                  f"latest on MAST: S{latest}; expected next (tess-point): {nxt or 'none listed'}")
    if w["period_err"] and w["t0_err"] and new:
        n = max(abs(r["epoch"]) for r in new)
        st = np.hypot(w["t0_err"], n * w["period_err"]) * 1440
        if st > 0.25 * w["duration_h"] * 60:
            report.append(f"- Note: the predicted times are uncertain by about {st:.0f} min at epoch {n}, a "
                          f"good fraction of the {w['duration_h'] * 60:.0f}-min transit; a low depth may be timing.")
    report.append("")
    report.append("| sector | epoch | product | depth ppm | error | SNR | |")
    report.append("|---|---|---|---|---|---|---|")
    for r in rows:
        report.append(f"| {r['sector']} | {r['epoch']} | {r['product']} | {r['depth_ppm']:.0f} | {r['err_ppm']:.0f} "
                      f"| {r['depth_ppm'] / r['err_ppm']:+.1f} | {'baseline' if r['baseline'] else 'new'} |")
    report.append("")
    if events:
        plot_events(star, w, events, plots)
        report.append(f"![{star}]({plots.name}/{slug(star)}.png)\n")
    logging.info(f"{star}: {v} - {txt}")


def plot_events(star, w, events, folder):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder.mkdir(parents=True, exist_ok=True)
    n = len(events)
    fig, axes = plt.subplots(n, 1, figsize=(8, 1.6 * n + 0.6), sharex=True, squeeze=False)
    t14 = w["duration_h"] / 24
    for ax, x in zip(axes[:, 0], events):
        dt = x["t"] - x["tc"]
        sel = np.abs(dt) < max(0.75, WINDOW_T14 * t14)
        h = dt[sel] * 24
        yy = x["y"][sel]
        base = np.polyval(np.polyfit(h[np.abs(h) > 0.75 * w["duration_h"]],
                                     yy[np.abs(h) > 0.75 * w["duration_h"]], 2), h)
        r = (yy / base - 1) * 1e6
        b = np.arange(h.min(), h.max(), 1 / 3)
        bm = [np.median(r[(h >= a) & (h < a + 1 / 3)]) if np.any((h >= a) & (h < a + 1 / 3)) else np.nan for a in b]
        ax.plot(h, r, ".", ms=1, color="0.75")
        ax.plot(b + 1 / 6, bm, "o-", ms=2.5, lw=1, color="C0")
        ax.plot(h, -w["depth_ppm"] * shape(h / 24, t14, TAU_FRAC * t14), color="C3", lw=1.2)
        ax.set_ylim(-4 * w["depth_ppm"], 3 * w["depth_ppm"])
        ax.text(0.01, 0.8, f"S{x['sector']} E{x['epoch']} {'(baseline)' if x['baseline'] else '(new)'}: "
                f"{x['depth_ppm']:.0f} +- {x['err_ppm']:.0f} ppm", transform=ax.transAxes, fontsize=8)
    axes[-1, 0].set_xlabel("hours from predicted mid-transit")
    fig.suptitle(f"{star}: transits measured this check (red: predicted dip)", fontsize=10)
    fig.tight_layout()
    fig.savefig(folder / f"{slug(star)}.png", dpi=90)
    plt.close(fig)


def check(con, stars):
    q = "SELECT * FROM watchlist" + (f" WHERE star IN ({','.join('?' * len(stars))})" if stars else "")
    con.row_factory = __import__("sqlite3").Row
    watched = [dict(r) for r in con.execute(q, stars)]
    con.row_factory = None
    if not watched:
        sys.exit("nothing to check (watchlist.py add ...)")
    day = datetime.date.today().isoformat()
    REPORTS.mkdir(parents=True, exist_ok=True)
    plots = REPORTS / day
    report = [f"# Watchlist check {day}\n"]
    for w in watched:
        try:
            check_star(con, w, report, plots)
        except Exception as e:
            logging.error(f"{w['star']}: check failed ({type(e).__name__}: {e}); rerun to retry")
            report.append(f"## {w['star']}: check failed ({type(e).__name__}: {e})\n")
    path = REPORTS / (f"{day}.md" if not stars else f"{day}_{'_'.join(slug(x) for x in stars)}.md")
    path.write_text("\n".join(report) + "\n")
    print(f"report: {path}")
    show(con)


def show(con):
    rows = con.execute("SELECT star, period, depth_ppm, verdict, summary, next_sectors, last_checked "
                       "FROM watchlist ORDER BY star").fetchall()
    for r in rows:
        print(f"{r[0]:15s} P={r[1]:.4f} d {r[2]:5.0f} ppm  {r[3] or 'not checked yet':13s} "
              f"{r[4] or ''}" + (f"  [next: {r[5]}]" if r[5] else "")
              + (f"  (checked {r[6][:10]})" if r[6] else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("add")
    a.add_argument("star")
    a.add_argument("--run"); a.add_argument("--candidate", type=int, default=1)
    a.add_argument("--period", type=float); a.add_argument("--t0", type=float)
    a.add_argument("--duration-h", type=float); a.add_argument("--depth-ppm", type=float)
    a.add_argument("--period-err", type=float); a.add_argument("--t0-err", type=float)
    a.add_argument("--sectors", help="baseline sectors (the ones the candidate was found in), comma-separated")
    a.add_argument("--reason")
    c = sub.add_parser("check")
    c.add_argument("stars", nargs="*")
    sub.add_parser("list")
    r = sub.add_parser("remove")
    r.add_argument("star")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    con = tessdb.connect()
    if args.cmd == "add":
        add(con, args.star, args)
    elif args.cmd == "check":
        check(con, args.stars)
    elif args.cmd == "list":
        show(con)
    elif args.cmd == "remove":
        with con:
            n = con.execute("DELETE FROM watchlist WHERE star = ?", (args.star,)).rowcount
        print(f"removed {args.star}" if n else f"{args.star} is not on the watchlist")


if __name__ == "__main__":
    main()
