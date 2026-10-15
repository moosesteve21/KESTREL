#!/usr/bin/env python3
"""
reduction_check.py - does the dip appear in every independent reduction of the data?

Usage:
    python reduction_check.py results/<name>/enriched.csv             # default selection
    python reduction_check.py results/<name>/enriched.csv --target "TIC 123456789" [--candidate 1]

A real transit is in the star's light, so it appears in every reduction of the
same observations. A processing artefact (for example from a systematics
correction near a momentum dump) appears in only one. For each selected
candidate this script loads every light curve MAST holds for the star:
  SPOC (2-min) and TESS-SPOC (FFI): PDCSAP (corrected) and SAP (raw) flux
  QLP: SAP (raw) and its own detrended flux (KSPSAP or DET)
and measures each predicted transit in each one: depth = mean flux just outside
minus inside, against a local straight-line baseline.

The reduction the search used (from the candidate's `sources`) is the reference.
For every other reduction, the events both cover are stacked and compared:
  supports     stacked depth > 3 sigma and >= half the reference depth
  contradicts  stacked depth < half the reference depth, by > 3 sigma
The candidate is flagged reduction_dependent when at least one independent
reduction contradicts it and none supports it.

Default selection: as pixel_check.py (clean, or double_period_variable with
enough transits, no nearby_variable_period_match). Writes
<out>/reduction_check.csv (each event in each reduction) and
<out>/reduction_summary.csv (one verdict per candidate).
"""
import argparse
import logging
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import lightkurve as lk                       # noqa: E402
from mastutil import retry, set_mast_timeout  # noqa: E402
set_mast_timeout()                    # give up on a stuck MAST query after 90 s, not 600

import pixel_check                            # default selection, shared with pixel_check.py

AUTHORS = ("SPOC", "TESS-SPOC", "QLP")
MIN_EXPTIME = 60                              # seconds; skip 20-s data
CORE = 0.4                                    # in-transit = |dt| < CORE * duration
NEAR = (0.75, 3.0)                            # baseline = this many durations either side...
MIN_BASE_DAYS = 0.5                           # ...but at least this far out
SUPPORT_SIGMA = 3.0
HALF = 0.5                                    # fraction of the reference depth
BKG_SIGMA, BKG_FRAC = 3.0, 0.3                # background_feature: the sky background rises in
                                              # transit by this many sigma and this fraction of the
                                              # depth (as a fraction of the star's flux): an
                                              # over-subtracted background event, not a transit
                                              # (a check PHT runs, Eisner et al. 2021)


def flux_columns(lc, author):
    """(label, column) pairs to measure for a product."""
    if author == "QLP":
        cols = [("QLP SAP", "sap_flux")]
        for c in ("kspsap_flux", "det_flux"):
            if c in lc.colnames:
                cols.append(("QLP detrended", c))
                break
        return cols
    return [(f"{author} PDCSAP", "pdcsap_flux"), (f"{author} SAP", "sap_flux")]


def primary_products(sources):
    """{sector: author} for the products the search used, from its sources string,
    e.g. 'SPOC/120s: 42,43,44; TESS-SPOC/200s: 70,71'."""
    out = {}
    for part in (sources or "").split(";"):
        m = re.match(r"\s*([\w-]+)/\d+s:\s*([\d,\s]+)", part)
        if m:
            for s in m.group(2).split(","):
                if s.strip():
                    out[int(s)] = m.group(1)
    return out


def event_depths(lc, col, P, T0, D):
    """Depth and error of each predicted transit in one flux column."""
    t = np.asarray(lc.time.value, float)
    v = np.asarray(lc[col].value, float)
    ok = np.isfinite(t) & np.isfinite(v)
    t, v = t[ok], v[ok]
    if len(t) < 50:
        return []
    out = []
    lo, hi = NEAR[0] * D, max(NEAR[1] * D, MIN_BASE_DAYS)
    for k in np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1):
        dt = t - (T0 + k * P)
        inn = np.abs(dt) < CORE * D
        base = (np.abs(dt) > lo) & (np.abs(dt) < hi)
        if inn.sum() < 3 or (base & (dt < 0)).sum() < 5 or (base & (dt > 0)).sum() < 5:
            continue                                        # need data on both sides
        c = np.polyfit(dt[base], v[base], 1)
        r = v / np.polyval(c, dt) - 1
        d = -np.mean(r[inn])
        e = np.std(r[base]) * np.sqrt(1 / inn.sum() + 1 / base.sum())
        out.append((int(k), float(d), float(e)))
    return out


def background_events(lc, P, T0, D):
    """How much the SAP background (SPOC, TESS-SPOC) rises during each predicted transit,
    as a fraction of the star's median SAP flux: (epoch, rise, error)."""
    if "sap_bkg" not in lc.colnames or "sap_flux" not in lc.colnames:
        return []
    t = np.asarray(lc.time.value, float)
    b = np.asarray(lc["sap_bkg"].value, float)
    f = np.asarray(lc["sap_flux"].value, float)
    ok = np.isfinite(t) & np.isfinite(b) & np.isfinite(f)
    t, b, f = t[ok], b[ok], f[ok]
    if len(t) < 50 or np.median(f) <= 0:
        return []
    scale = np.median(f)
    out = []
    lo, hi = NEAR[0] * D, max(NEAR[1] * D, MIN_BASE_DAYS)
    for k in np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1):
        dt = t - (T0 + k * P)
        inn = np.abs(dt) < CORE * D
        base = (np.abs(dt) > lo) & (np.abs(dt) < hi)
        if inn.sum() < 3 or (base & (dt < 0)).sum() < 5 or (base & (dt > 0)).sum() < 5:
            continue
        c = np.polyfit(dt[base], b[base], 1)
        r = (b - np.polyval(c, dt)) / scale
        out.append((int(k), float(np.mean(r[inn])),
                    float(np.std(r[base]) * np.sqrt(1 / inn.sum() + 1 / base.sum()))))
    return out


def stack(ev):
    """Inverse-variance mean of (depth, error) pairs."""
    d = np.array([x[0] for x in ev])
    e = np.array([x[1] for x in ev])
    w = 1 / e ** 2
    return float(np.sum(d * w) / w.sum()), float(1 / np.sqrt(w.sum()))


def check(r):
    name, k = r["target"], int(r["candidate"])
    P, T0, D = float(r["period_d"]), float(r["t0_btjd"]), float(r["duration_h"]) / 24
    prim = primary_products(r.get("sources", ""))
    logging.info(f"{name} #{k}: P = {P:.5f} d, dur = {D * 24:.2f} h; search used "
                 f"{', '.join(f'S{s} {a}' for s, a in sorted(prim.items())) or 'unknown'}")
    sr = retry(lambda: lk.search_lightcurve(name, author=AUTHORS), f"{name}: MAST search")
    rows, bkg = [], []
    for item in sr:
        author = str(item.table["author"][0])
        sector = int(item.table["sequence_number"][0])
        if float(item.table["exptime"][0]) < MIN_EXPTIME:
            continue
        try:
            lc = item.download(quality_bitmask="default")
        except Exception as e:
            logging.warning(f"  {author} S{sector}: download failed ({type(e).__name__})")
            continue
        if author == "QLP" and "quality" in lc.colnames:   # QLP's own flag bits
            q = np.asarray(lc.quality.value if hasattr(lc.quality, "value") else lc.quality)
            lc = lc[q.astype(np.int64) == 0]
        if author in ("SPOC", "TESS-SPOC"):
            bkg += background_events(lc, P, T0, D)
        for label, col in flux_columns(lc, author):
            if col not in lc.colnames:
                continue
            is_primary = prim.get(sector) == author and col == ("sap_flux" if author == "QLP"
                                                                 else "pdcsap_flux")
            for epoch, d, e in event_depths(lc, col, P, T0, D):
                rows.append(dict(target=name, candidate=k, sector=sector, reduction=label,
                                 primary=is_primary, epoch=epoch, depth_ppm=d * 1e6, err_ppm=e * 1e6))
    if not rows:
        logging.info("  -> no reduction covers a transit")
        return rows, dict(target=name, candidate=k, verdict="no data", detail="")
    df = pd.DataFrame(rows)
    ref = df[df.primary]
    if ref.empty:
        logging.info("  -> the search's own reduction isn't available; nothing to compare")
        return rows, dict(target=name, candidate=k, verdict="no reference", detail="")
    ref_d, ref_e = stack(list(zip(ref.depth_ppm, ref.err_ppm)))
    logging.info(f"  reference (search's reduction): {len(ref)} events, stacked {ref_d:.0f} +- {ref_e:.0f} ppm")
    n_sup = n_con = 0
    parts = []
    others = df[~df.primary]
    for label in sorted(others.reduction.unique()):
        o = others[others.reduction == label]
        common = sorted(set(zip(o.sector, o.epoch)) & set(zip(ref.sector, ref.epoch)))
        if not common:
            continue
        oc = o[[(s, ep) in common for s, ep in zip(o.sector, o.epoch)]]
        rc = ref[[(s, ep) in common for s, ep in zip(ref.sector, ref.epoch)]]
        d, e = stack(list(zip(oc.depth_ppm, oc.err_ppm)))
        rd, _ = stack(list(zip(rc.depth_ppm, rc.err_ppm)))
        if d > SUPPORT_SIGMA * e and d >= HALF * rd:
            word = "supports"
            n_sup += 1
        elif d < HALF * rd - SUPPORT_SIGMA * e:
            word = "contradicts"
            n_con += 1
        else:
            word = "inconclusive"
        parts.append(f"{label}: {d:.0f}+-{e:.0f} ppm vs {rd:.0f} on {len(common)} events, {word}")
        logging.info(f"  {label:18s} {d:6.0f} +- {e:4.0f} ppm (reference {rd:6.0f} on the same "
                     f"{len(common)} events): {word}")
    if n_con >= 1 and n_sup == 0:
        verdict = "reduction_dependent"
    elif n_sup >= 1 and n_con == 0:
        verdict = "consistent"
    elif n_sup == 0 and n_con == 0:
        verdict = "inconclusive"
    else:
        verdict = "mixed"
    bkg_d = bkg_e = None
    if bkg:
        bd, be = stack([(d, e) for _, d, e in bkg if e > 0])
        bkg_d, bkg_e = round(bd * 1e6), round(be * 1e6)
        feature = bd * 1e6 > BKG_SIGMA * be * 1e6 and bd * 1e6 > BKG_FRAC * ref_d
        parts.append(f"background rise in transit {bkg_d}+-{bkg_e} ppm of the star's flux over {len(bkg)} events"
                     + ("; background_feature" if feature else ""))
        logging.info(f"  background in transit: {bkg_d:+d} +- {bkg_e} ppm of the star's flux"
                     + (" -> background_feature" if feature else ""))
    logging.info(f"  -> {verdict} ({n_sup} supporting, {n_con} contradicting)")
    return rows, dict(target=name, candidate=k, verdict=verdict, n_support=n_sup,
                      n_contradict=n_con, reference_ppm=round(ref_d), bkg_rise_ppm=bkg_d,
                      bkg_err_ppm=bkg_e, detail="; ".join(parts))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", help="enriched.csv or candidates_all.csv")
    ap.add_argument("--target", help="check this star only")
    ap.add_argument("--candidate", type=int, help="with --target: this candidate only")
    ap.add_argument("--out", help="output folder (default: reduction_check beside the input)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("lightkurve").setLevel(logging.WARNING)
    src = Path(args.candidates)
    outdir = Path(args.out) if args.out else src.parent / "reduction_check"
    outdir.mkdir(parents=True, exist_ok=True)
    df = pixel_check.select_rows(pd.read_csv(src), args.target, args.candidate)
    logging.info(f"{len(df)} candidate(s) to check")
    rows, summary = [], []
    for _, r in df.iterrows():
        try:
            ev, verdict = check(r)
            rows += ev
            summary.append(verdict)
        except Exception:
            logging.exception(f"  {r['target']} #{r['candidate']} failed")
    if rows:
        pd.DataFrame(rows).to_csv(outdir / "reduction_check.csv", index=False)
    if summary:
        pd.DataFrame(summary).to_csv(outdir / "reduction_summary.csv", index=False)
        logging.info(f"wrote {outdir / 'reduction_summary.csv'}")


if __name__ == "__main__":
    main()
