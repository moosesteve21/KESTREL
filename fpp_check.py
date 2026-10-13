#!/usr/bin/env python3
"""
fpp_check.py - TRICERATOPS false-positive probability for candidates whose transit
fit is consistent with their star.

Usage:
    python fpp_check.py results/<name>/enriched.csv
    python fpp_check.py results/<name>/enriched.csv --target "TIC 302453503" [--candidate 1]
    options: --draws N (default 1000000), --fits DIR (default: transit_fit beside the
             input), --out DIR (default: fpp_check beside the input)

Runs after transit_fit.py, on the candidates it judged consistent, using the same
normalised transits. Each runs vetting/fpp_triceratops.py in the vetting
environment (.venv-vet), which queries Gaia and MAST: retried up to 3 times.

Label (thresholds of Giacalone et al. 2021):
  validation_level   FPP < 0.015 and NFPP < 0.001
  likely_planet      FPP < 0.5 and NFPP < 0.001
  likely_nearby_fp   NFPP > 0.1 (another star is the likely source)
  unclear            anything else
  no_aperture        no SPOC or TESS-SPOC sector (QLP publishes no apertures)
It records the result, it doesn't decide: no contrast curve is used, so the
unresolved-companion scenarios stay open until high-resolution imaging exists.
Writes <out>/fpp.csv.
"""
import argparse
import json
import logging
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
VET_PY = HERE / ".venv-vet" / "bin" / "python"
TRIES = 3


def sectors_of(sources, author):
    """Sectors of one author in a sources string like 'SPOC/120s: 1,28; TESS-SPOC/200s: 68'."""
    out = []
    for part in str(sources).split(";"):
        m = re.match(r"\s*([\w-]+)/\d+s:\s*([\d,\s]+)", part)
        if m and m.group(1) == author:
            out += [int(x) for x in m.group(2).split(",") if x.strip()]
    return out


def label(fpp, nfpp):
    if fpp < 0.015 and nfpp < 0.001:
        return "validation_level"
    if fpp < 0.5 and nfpp < 0.001:
        return "likely_planet"
    if nfpp > 0.1:
        return "likely_nearby_fp"
    return "unclear"


def check_one(fit, outdir, fitdir, draws):
    name, k = fit["target"], int(fit["candidate"])
    base = dict(target=name, candidate=k)
    m = re.fullmatch(r"TIC\s*(\d+)", name)
    if not m:
        return dict(base, label="error", detail="needs a TIC id")
    par = dict(tic=int(m.group(1)), P=fit["period"], T0=fit["t0"], depth=fit["depth_ppm"] * 1e-6,
               t14_h=fit["t14_h"], events=str((fitdir / fit["events_file"]).resolve()), draws=draws,
               spoc_sectors=sectors_of(fit["sources"], "SPOC"),
               tess_spoc_sectors=sectors_of(fit["sources"], "TESS-SPOC"),
               workdir=str((outdir / "work").resolve()))
    if not par["spoc_sectors"] and not par["tess_spoc_sectors"]:
        return dict(base, label="no_aperture", detail="QLP-only data")
    pj = outdir / f"{re.sub(r'[^\w.-]+', '_', name)}_c{k}_fpp.json"
    pj.write_text(json.dumps(par))
    for attempt in range(1, TRIES + 1):
        logging.info(f"{name} #{k}: TRICERATOPS, {draws} draws (attempt {attempt})")
        res = subprocess.run([str(VET_PY), str(HERE / "vetting" / "fpp_triceratops.py"), str(pj.resolve())],
                             capture_output=True, text=True)
        lines = [l for l in res.stdout.splitlines() if l.startswith("{")]
        if res.returncode == 0 and lines:
            r = json.loads(lines[-1])
            if "error" in r:
                return dict(base, label="no_aperture", detail=r["error"])
            lab = label(r["fpp"], r["nfpp"])
            logging.info(f"  -> FPP {r['fpp']:.4f}, NFPP {r['nfpp']:.4f}: {lab}")
            return dict(base, label=lab, fpp=round(r["fpp"], 5), nfpp=round(r["nfpp"], 5),
                        prob_tp=round(r["prob_tp"], 4), draws=draws,
                        sectors=",".join(map(str, r["sectors"])),
                        scenarios="; ".join(f"{k_}={v}" for k_, v in r["scenarios"].items()))
        err = (res.stderr.strip().splitlines() or ["?"])[-1][:200]
        logging.warning(f"  attempt {attempt} failed: {err}")
        if attempt < TRIES:
            time.sleep(60)
    return dict(base, label="error", detail=err)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", help="enriched.csv (used for its folder and the selection)")
    ap.add_argument("--target")
    ap.add_argument("--candidate", type=int)
    ap.add_argument("--draws", type=int, default=1_000_000)
    ap.add_argument("--fits")
    ap.add_argument("--out")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    src = Path(args.candidates)
    fitdir = Path(args.fits) if args.fits else src.parent / "transit_fit"
    outdir = Path(args.out) if args.out else src.parent / "fpp_check"
    outdir.mkdir(parents=True, exist_ok=True)
    if not (fitdir / "transit_fit.csv").exists():
        logging.info(f"no {fitdir / 'transit_fit.csv'}: nothing to check")
        return 0
    fits_df = pd.read_csv(fitdir / "transit_fit.csv")
    sel = fits_df[fits_df.verdict == "consistent"]
    if args.target:
        sel = sel[sel.target == args.target]
        if args.candidate:
            sel = sel[sel.candidate == args.candidate]
    logging.info(f"{len(sel)} candidate(s) with a consistent transit fit")
    if sel.empty:
        return 0
    rows = []
    for _, f in sel.iterrows():
        try:
            rows.append(check_one(f, outdir, fitdir, args.draws))
        except Exception as e:
            logging.exception(f"  {f['target']} #{f['candidate']} failed")
            rows.append(dict(target=f["target"], candidate=int(f["candidate"]), label="error",
                             detail=f"{type(e).__name__}: {e}"[:200]))
    pd.DataFrame(rows).to_csv(outdir / "fpp.csv", index=False)
    logging.info(f"wrote {outdir / 'fpp.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
