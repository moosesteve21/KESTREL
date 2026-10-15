#!/usr/bin/env python3
"""
depth_map_check.py - per-source pixel depth map for clean candidates (depth_map.py).

Usage:
    python depth_map_check.py results/<name>/enriched.csv
    python depth_map_check.py results/<name>/enriched.csv --target "TIC 123456789" [--candidate 1]
    options: --out DIR (default: depth_map beside the input), --fits DIR (default:
             transit_fit beside the input)

Runs after transit_fit.py: the ephemeris is the transit fit's when it is consistent
(over a long baseline the search's period can drift by tens of minutes), else the
search's, with the shape fit's duration. Every clean candidate whose name is a TIC.

Verdict, from the chi-square of each source's star image against the deficit map,
summed over sectors:
  on_target      the target fits within TARGET_DCHI2 of the best source. Others that
                 fit as well are listed in detail (the map can't split close stars)
  off_target     another star fits better by more than OFF_DCHI2 and is the best source
                 in at least half the sectors (flag: depth_map_off_target)
  inconclusive   neither; or the verdict flips between on and off target when any one
                 sector is left out (detail names it: e.g. a sector dominated by
                 scattered light, which DEPTHMAP_SKIP leaves out by hand)
  not_detected   no source beats "no dip anywhere" by DETECT_DCHI2 (too shallow to
                 locate in the FFI cutouts: not evidence against the target)
  no_data / error
fpp_check.py skips off_target candidates. Writes <out>/depth_map.csv and, per
candidate, <star>_c<k>.png and .txt.
"""
import argparse
import logging
import re
import sys
import time
from pathlib import Path

import pandas as pd

import depth_map as dm

TARGET_DCHI2 = 9.0          # target within this of the best source: on target
OFF_DCHI2 = 25.0            # another star better by more than this (and wins half the sectors): off target
DETECT_DCHI2 = 16.0         # the best source must beat "no dip" by this (~4 sigma)
TRIES = 3


def ephemeris(r, fits):
    """(P, T0, T14 in days, source) for a candidate row."""
    dur = r.get("dur_fit_h")
    t14 = float(dur if pd.notna(dur) and str(dur) != "" else r["duration_h"]) / 24
    if fits is not None:
        f = fits[(fits.target == r["target"]) & (fits.candidate == int(r["candidate"]))]
        if len(f) and f.iloc[0]["verdict"] == "consistent" and pd.notna(f.iloc[0].get("period")):
            f = f.iloc[0]
            return float(f["period"]), float(f["t0"]), float(f["t14_h"]) / 24, "transit fit"
        if len(f) and f.iloc[0]["verdict"] in ("density_mismatch", "grazing_density_uncertain") \
                and pd.notna(f.iloc[0].get("period")):
            f = f.iloc[0]
            return float(f["period"]), float(f["t0"]), t14, "transit fit (search duration)"
    return float(r["period_d"]), float(r["t0_btjd"]), t14, "search"


def verdict(res):
    """(verdict, source, target dchi2, detection dchi2, sectors won, detail), checked
    by leaving each sector out in turn."""
    v = placement(res)
    secs = res["sectors"]
    if v[0] not in ("on_target", "off_target", "inconclusive") or len(secs) < 2:
        return v
    flips = []
    for i, s in enumerate(secs):
        rest = secs[:i] + secs[i + 1:]
        chi = {k: x - s[5].get(k, s[6]) for k, x in res["chi"].items()}
        w = placement(dict(res, chi=chi, chi_null=res["chi_null"] - s[6], sectors=rest))
        if w[0] != v[0] and {w[0], v[0]} & {"on_target", "off_target"} and \
                {w[0], v[0]} <= {"on_target", "off_target", "inconclusive"}:
            flips.append(f"without S{s[0]}: {w[0]}" + (f" ({w[1]})" if w[0] == "off_target" else ""))
    if not flips:
        return v
    if v[0] == "inconclusive":
        return v[:5] + ("; ".join(x for x in [v[5]] + flips if x),)
    return ("inconclusive", v[1], v[2], v[3], v[4],
            f"{v[0]} with every sector, but " + "; ".join(flips))


def placement(res):
    """Verdict from the summed chi-squares."""
    tic, chi, secs = str(res["tic"]), res["chi"], res["sectors"]
    if not secs:
        return "no_data", "", None, None, "", "no sector with a covered transit"
    best = min(chi, key=chi.get)
    t_d = chi[tic] - chi[best]
    det = res["chi_null"] - chi[best]
    sep = {s[0]: (s[3], s[4]) for s in res["stars"]}
    if det < DETECT_DCHI2:
        return "not_detected", "", round(t_d, 1), round(det, 1), "", "no source beats 'no dip' clearly"
    if t_d <= TARGET_DCHI2:
        alt = sorted((chi[s] - chi[best], s) for s in chi if s != tic and chi[s] - chi[best] <= TARGET_DCHI2)
        detail = ("also fits: " + "; ".join(f"TIC {s} (T={sep[s][0]:.1f}, {sep[s][1]:.0f}\")" for _, s in alt[:3])
                  if alt else "")
        return "on_target", f"TIC {tic}", round(t_d, 1), round(det, 1), "", detail
    wins = sum(1 for s in secs if s[2] == best)
    won = f"{wins}/{len(secs)}"
    src = f"TIC {best} (T={sep[best][0]:.1f}, {sep[best][1]:.0f}\")"
    if t_d > OFF_DCHI2 and wins >= len(secs) / 2:
        return "off_target", src, round(t_d, 1), round(det, 1), won, f"best source in {won} sectors"
    return "inconclusive", src, round(t_d, 1), round(det, 1), won, f"best source in {won} sectors"


def check_one(r, fits, outdir):
    name, k = r["target"], int(r["candidate"])
    base = dict(target=name, candidate=k)
    m = re.fullmatch(r"TIC\s*(\d+)", str(name))
    if not m:
        return dict(base, verdict="error", detail="needs a TIC id")
    P, T0, T14, src = ephemeris(r, fits)
    stem = outdir / f"{re.sub(r'[^\w.-]+', '_', name)}_c{k}"
    for attempt in range(1, TRIES + 1):
        lines = []
        try:
            logging.info(f"{name} #{k}: P {P:.6f}, T0 {T0:.5f}, T14 {T14 * 24:.2f} h ({src})")
            res = dm.depth_map(int(m[1]), P, T0, T14, say=lines.append)
            break
        except Exception as e:
            logging.warning(f"  attempt {attempt} failed: {type(e).__name__}: {e}")
            if attempt == TRIES:
                return dict(base, verdict="error", detail=f"{type(e).__name__}: {e}"[:200])
            time.sleep(60)
    lines.append("")
    dm.ranking(res, say=lines.append)
    stem.with_suffix(".txt").write_text(f"ephemeris ({src}): P {P:.6f} T0 {T0:.5f} T14 {T14 * 24:.2f} h\n"
                                        + "\n".join(lines) + "\n")
    dm.plot(res, stem.with_suffix(".png"))
    v, source, t_d, det, won, detail = verdict(res)
    logging.info(f"  -> {v}: {source} (target dchi2 {t_d}, detection {det}) {detail}")
    return dict(base, verdict=v, source=source, target_dchi2=t_d, detect_dchi2=det, sectors_won=won,
                n_sectors=len(res["sectors"]), ephemeris=src, detail=detail)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", help="enriched.csv")
    ap.add_argument("--target")
    ap.add_argument("--candidate", type=int)
    ap.add_argument("--out")
    ap.add_argument("--fits")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    src = Path(args.candidates)
    outdir = Path(args.out) if args.out else src.parent / "depth_map"
    fitcsv = (Path(args.fits) if args.fits else src.parent / "transit_fit") / "transit_fit.csv"
    fits = pd.read_csv(fitcsv) if fitcsv.exists() else None
    df = pd.read_csv(src)
    if args.target:
        sel = df[df.target == args.target]
        if args.candidate:
            sel = sel[sel.candidate == args.candidate]
    else:
        sel = df[df.status == "clean"]
    logging.info(f"{len(sel)} candidate(s) for a depth map")
    if sel.empty:
        if not args.target:
            (outdir / "depth_map.csv").unlink(missing_ok=True)   # don't leave an earlier run's results
        return 0
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, r in sel.iterrows():
        try:
            rows.append(check_one(r, fits, outdir))
        except Exception as e:
            logging.exception(f"  {r['target']} #{r['candidate']} failed")
            rows.append(dict(target=r["target"], candidate=int(r["candidate"]), verdict="error",
                             detail=f"{type(e).__name__}: {e}"[:200]))
    pd.DataFrame(rows).to_csv(outdir / "depth_map.csv", index=False)
    logging.info(f"wrote {outdir / 'depth_map.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
