"""
investigations/other_pipelines.py - is a candidate's dip in light curves from
pipelines the search doesn't use?

Usage:
    python investigations/other_pipelines.py TIC P_DAYS T0_BTJD DURATION_H [OUTDIR]

Every TESS light curve MAST holds for the star, except the search's own (SPOC,
TESS-SPOC, QLP): lightkurve-readable ones (GSFC-ELEANOR-LITE, TGLC, ...) through
lightkurve, and TARS (time, flux only) and TASOC (FLUX_CORR, relative ppm) read
directly. Each predicted transit's depth is measured as in reduction_check.py
(mean outside minus inside, against a local straight line), per sector and
stacked per pipeline. The stack's significance is also judged against stacks of
the same size at random phases, which includes red noise.
"""
import sys
import warnings
from glob import glob
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import lightkurve as lk
import numpy as np
from astropy.io import fits
from astroquery.mast import Observations

from mastutil import retry, set_mast_timeout
from reduction_check import event_depths, stack

set_mast_timeout()
SEARCH_AUTHORS = {"SPOC", "TESS-SPOC", "QLP"}


def report(label, per_sector, rand, rng):
    ev = [x for v in per_sector.values() for x in v]
    for sec, v in sorted(per_sector.items()):
        if v:
            d, e = stack(v)
            print(f"  {label} S{sec}: {len(v)} transits, {d * 1e6:.0f} +- {e * 1e6:.0f} ppm")
    if not ev:
        print(f"{label}: no transits covered")
        return
    d, e = stack(ev)
    rd = np.array([x[0] for x in rand])
    emp = np.std([np.mean(rng.choice(rd, len(ev))) for _ in range(2000)]) if len(rd) >= len(ev) else np.nan
    print(f"{label} (all): {len(ev)} transits, stacked {d * 1e6:.0f} +- {e * 1e6:.0f} ppm "
          f"({d / e:.1f} sigma); against random phases {d / emp:.1f} sigma")


def measure(lc, P, T0, D, rng, per_sector, rand, sec):
    e1 = event_depths(lc, "flux", P, T0, D)
    per_sector.setdefault(sec, []).extend((d, e) for _, d, e in e1)
    for off in rng.uniform(0.15, 0.85, 8):
        rand.extend((d, e) for _, d, e in event_depths(lc, "flux", P, T0 + off * P, D))


def main():
    tic, P, T0, D = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) / 24
    out = Path(sys.argv[5]) if len(sys.argv) > 5 else Path("scratch/hlsp") / tic
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    name = f"TIC {tic}"
    sr = retry(lambda: lk.search_lightcurve(name), f"{name}: MAST search")
    authors = sorted(set(sr.table["author"]) - SEARCH_AUTHORS) if len(sr) else []
    print(f"{name}: pipelines besides the search's own: {', '.join(authors) or 'none'}")
    direct = {"TARS", "TASOC"}
    for author in [a for a in authors if a not in direct]:
        per_sector, rand = {}, []
        for item in sr[sr.table["author"] == author]:
            sec = int(item.table["sequence_number"][0])
            try:
                lc = item.download().normalize()
            except Exception as e:
                print(f"  {author} S{sec}: not readable ({type(e).__name__})")
                continue
            measure(lc, P, T0, D, rng, per_sector, rand, sec)
        report(author, per_sector, rand, rng)
    # TARS and TASOC: download and read the files directly
    want = [a for a in authors if a in direct]
    if want:
        obs = retry(lambda: Observations.query_criteria(target_name=str(int(tic)), provenance_name=want),
                    f"{name}: HLSP query")
        prods = Observations.get_product_list(obs)
        keep = [str(f).endswith("_lc.fits") and "tars" in str(f) or str(f).endswith("cbv-lc.fits")
                for f in prods["productFilename"]]
        if any(keep):
            Observations.download_products(prods[keep], download_dir=str(out))
        for label, pattern, col in (("TARS", "**/hlsp_tars_*_lc.fits", "FLUX"),
                                    ("TASOC", "**/hlsp_tasoc_*cbv-lc.fits", "FLUX_CORR")):
            if label not in want:
                continue
            per_sector, rand = {}, []
            for path in sorted(glob(str(out / pattern), recursive=True)):
                with fits.open(path) as h:
                    d = h[1].data
                    t, f = np.asarray(d["TIME"], float), np.asarray(d[col], float)
                    q = np.asarray(d["QUALITY"]) if "QUALITY" in d.names else np.zeros(len(t))
                    sec = h[0].header.get("SECTOR")
                ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
                t, f = t[ok], f[ok]
                if len(t) < 50:
                    continue
                f = 1 + f / 1e6 if abs(np.median(f)) < 0.1 * np.std(f) + 1e-9 else f / np.median(f)
                measure(lk.LightCurve(time=t, flux=f), P, T0, D, rng, per_sector, rand, sec)
            report(label, per_sector, rand, rng)


if __name__ == "__main__":
    main()
