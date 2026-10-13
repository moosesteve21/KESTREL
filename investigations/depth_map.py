#!/usr/bin/env python3
"""Where on the sky is a candidate's dip? A per-pixel depth map from TESScut cutouts,
for stars without 2-min data (no SPOC centroids). Generalised from tic240488148/depth_map.py.

For each TESScut sector: every pixel's light curve gets a local linear baseline around
each transit, and the in-transit flux deficit is averaged over all transits, giving a
deficit map (e-/s) with errors from the transit-to-transit scatter. Each candidate
source (the target and its brighter neighbours) predicts that map as a star image at
its position: a 2-D Gaussian whose width is fitted to the target's own out-of-transit
image. A free amplitude and offset are fitted per sector; the chi-square of each
hypothesis, summed over sectors, says which star the dip comes from.

Usage: .venv/bin/python investigations/depth_map.py TIC P_DAYS T0_BTJD T14_HOURS OUTDIR
Writes depth_map.png and depth_map.txt in OUTDIR. DEPTHMAP_SKIP=63,64 leaves sectors out
(for example one dominated by scattered light).
"""
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import astropy.units as u
import lightkurve as lk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.mast import Catalogs
from scipy.optimize import least_squares

import sys
TIC = int(sys.argv[1])
P, T0, T14 = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) / 24
OUT = Path(sys.argv[5]); OUT.mkdir(parents=True, exist_ok=True)
SIZE = 15
WIN = 3 * T14                      # baseline window each side of mid-transit
MIN_TRANSITS = 1                   # per sector
import os
SKIP = {int(x) for x in os.environ.get("DEPTHMAP_SKIP", "").split(",") if x.strip()}   # sectors to leave out
lines = []


def say(s=""):
    print(s)
    lines.append(s)


tgt = Catalogs.query_criteria(catalog="TIC", ID=TIC)[0]
coord = SkyCoord(float(tgt["ra"]) * u.deg, float(tgt["dec"]) * u.deg)
nb = Catalogs.query_region(coord, radius=2.5 * u.arcmin, catalog="TIC").to_pandas()
nb = nb[(nb["Tmag"] < float(tgt["Tmag"]) + 5.5)].sort_values("Tmag")
stars = [(str(r["ID"]), float(r["ra"]), float(r["dec"]), float(r["Tmag"]), float(r["dstArcSec"]))
         for _, r in nb.iterrows()]
say(f"TIC {TIC}: T = {float(tgt['Tmag']):.2f}; {len(stars) - 1} neighbours within 2.5' brighter than T = "
    f"{float(tgt['Tmag']) + 5.5:.1f}")


def gauss(x, y, x0, y0, s):
    return np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * s ** 2))


sr = lk.search_tesscut(f"TIC {TIC}")
maps, chi = [], {s[0]: 0.0 for s in stars}
chi_null = 0.0
yy, xx = np.mgrid[0:SIZE, 0:SIZE]
for row in sr:
    sector = int(row.table["sequence_number"][0])
    if sector in SKIP:
        say(f"S{sector}: left out (DEPTHMAP_SKIP)")
        continue
    try:
        tpf = row.download(cutout_size=SIZE, quality_bitmask="hard")
    except Exception as e:
        say(f"S{sector}: download failed ({type(e).__name__})")
        continue
    t = np.asarray(tpf.time.value)
    cube = np.asarray(tpf.flux.value, float)
    ok = np.isfinite(t) & np.all(np.isfinite(cube), axis=(1, 2))
    t, cube = t[ok], cube[ok]
    deficits, errs = [], []
    for k in np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1):
        tc = T0 + k * P
        inn = np.abs(t - tc) < 0.35 * T14
        oot = (np.abs(t - tc) > 0.75 * T14) & (np.abs(t - tc) < WIN)
        if inn.sum() < 1 or (oot & (t < tc)).sum() < 3 or (oot & (t > tc)).sum() < 3:
            continue
        A = np.column_stack([np.ones(oot.sum()), t[oot] - tc])
        coef, *_ = np.linalg.lstsq(A, cube[oot].reshape(oot.sum(), -1), rcond=None)
        base = np.column_stack([np.ones(inn.sum()), t[inn] - tc]) @ coef
        deficits.append((base - cube[inn].reshape(inn.sum(), -1)).mean(axis=0).reshape(SIZE, SIZE))
        # this transit's per-pixel error from the out-of-transit scatter about the baseline
        resid = cube[oot].reshape(oot.sum(), -1) - A @ coef
        errs.append((resid.std(axis=0, ddof=2) / np.sqrt(inn.sum())).reshape(SIZE, SIZE))
    if len(deficits) < MIN_TRANSITS:
        say(f"S{sector}: no transit covered; skipped")
        continue
    D = np.array(deficits)
    dmap = D.mean(axis=0)
    # errors: transit-to-transit scatter where there are enough transits, else (and as a
    # floor) the out-of-transit scatter propagated, so long periods with 1-2 transits work
    err_oot = np.sqrt(np.sum(np.array(errs) ** 2, axis=0)) / len(D)
    err = np.maximum(D.std(axis=0, ddof=1) / np.sqrt(len(D)), err_oot) if len(D) >= 4 else err_oot
    img = np.median(cube, axis=0)
    # star positions in this cutout
    pos = {}
    for sid, ra, dec, *_ in stars:
        x, y = tpf.wcs.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
        pos[sid] = (float(x), float(y))
    # PSF width and a position correction from the target's out-of-transit image
    x0, y0 = pos[str(TIC)]
    sel = gauss(xx, yy, x0, y0, 2.0) > 0.05

    def r_img(p):
        a, xc, yc, s, c = p
        return ((a * gauss(xx, yy, xc, yc, s) + c) - img)[sel]

    fit = least_squares(r_img, [img[sel].max(), x0, y0, 1.0, np.median(img)],
                        bounds=([0, x0 - 1, y0 - 1, 0.4, -np.inf], [np.inf, x0 + 1, y0 + 1, 3, np.inf]))
    _, xc, yc, sig, _ = fit.x
    dx, dy = xc - x0, yc - y0                       # shift every star by the same WCS offset
    # hypotheses: dip from star s -> deficit map = A * gauss at s + c
    use = np.hypot(xx - xc, yy - yc) < 5            # the region around the target
    w = 1 / err[use]
    c2 = {}
    for sid, (px, py) in pos.items():
        g = gauss(xx, yy, px + dx, py + dy, sig)[use]
        X = np.column_stack([g, np.ones_like(g)]) * w[:, None]
        coef, *_ = np.linalg.lstsq(X, dmap[use] * w, rcond=None)
        if coef[0] < 0:                             # a brightening isn't a transit
            coef[0] = 0
            coef[1] = np.sum(dmap[use] * w ** 2) / np.sum(w ** 2)
        c2[sid] = float(np.sum((dmap[use] * w - X @ coef) ** 2))
    c0 = float(np.sum(((dmap[use] - np.sum(dmap[use] * w ** 2) / np.sum(w ** 2)) * w) ** 2))
    for sid in chi:
        chi[sid] += c2.get(sid, c0)
    chi_null += c0
    best = min(c2, key=c2.get)
    say(f"S{sector}: {len(D)} transits, PSF sigma {sig:.2f} px; best source TIC {best} "
        f"(target dchi2 {c2[str(TIC)] - c2[best]:+.1f}; no-dip dchi2 {c0 - c2[str(TIC)]:+.1f} for the target)")
    maps.append((sector, dmap, err, img, {k: (v[0] + dx, v[1] + dy) for k, v in pos.items()}))

say()
say("Summed over sectors, chi-square relative to the best hypothesis (lower = better):")
best = min(chi, key=chi.get)
for sid, ra, dec, tmag, sep in stars:
    say(f"  TIC {sid:>10}  T={tmag:5.2f}  {sep:5.1f}\"   dchi2 = {chi[sid] - chi[best]:7.1f}"
        + ("   <- target" if sid == str(TIC) else ""))
say(f"  no dip anywhere                 dchi2 = {chi_null - chi[best]:7.1f}")
(OUT / "depth_map.txt").write_text("\n".join(lines) + "\n")

fig, axes = plt.subplots(2, len(maps), figsize=(3.2 * len(maps), 6.4), squeeze=False)
for i, (sector, dmap, err, img, pos) in enumerate(maps):
    for ax, data, title in ((axes[0, i], img, f"S{sector} image"),
                            (axes[1, i], (dmap - np.median(dmap)) / np.nanmedian(err), f"S{sector} deficit / error")):
        ax.imshow(data, origin="lower", cmap="viridis" if "image" in title else "RdBu_r",
                  vmin=None if "image" in title else -4, vmax=None if "image" in title else 4)
        for sid, (px, py) in pos.items():
            ax.plot(px, py, "x" if sid != str(TIC) else "*", color="w" if sid != str(TIC) else "gold",
                    ms=6 if sid != str(TIC) else 12)
        ax.set_title(title, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f"TIC {TIC}: in-transit flux deficit per pixel (gold star = target, x = neighbours)")
fig.tight_layout()
fig.savefig(OUT / "depth_map.png", dpi=90)
