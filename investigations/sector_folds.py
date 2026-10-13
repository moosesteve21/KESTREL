#!/usr/bin/env python3
"""A candidate's transit folded separately in each sector (for Talk posts).
Each transit gets a local linear baseline (out of transit, +-3 T14); points are binned
to 15 minutes.
Usage: .venv/bin/python investigations/sector_folds.py TIC P_DAYS T0_BTJD T14_HOURS DEPTH_PPM OUTDIR"""
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import lightkurve as lk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import find_planets as fp

TIC = sys.argv[1]
P, T0, T14, DEPTH = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) / 24, float(sys.argv[5]) * 1e-6
OUT = Path(sys.argv[6]); OUT.mkdir(parents=True, exist_ok=True)
sr = lk.search_lightcurve(f"TIC {TIC}", author=fp.AUTHORS)
rows, _ = fp.choose_products(sr)
lcs = [(int(sr.table["sequence_number"][i]), f"{sr.table['author'][i]} {int(sr.table['exptime'][i])} s",
        fp.clean(sr[i].download(quality_bitmask="default"))) for i in rows]
lcs = [x for x in lcs if x[2] is not None]
ncol = 2 if len(lcs) == 4 else min(3, len(lcs)); nrow = -(-len(lcs) // ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.3 * nrow), sharex=True, sharey=True, squeeze=False)
for ax, (sec, prod, lc) in zip(axes.flat, lcs):
    t, y = fp.arr(lc.time), fp.arr(lc.flux)
    hs, rs = [], []
    for k in np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1):
        tc = T0 + k * P
        w = np.abs(t - tc) < 3 * T14
        oot = w & (np.abs(t - tc) > 0.75 * T14)
        if oot.sum() < 6 or np.sum(np.abs(t - tc) < T14 / 2) < 1:
            continue
        cf = np.polyfit(t[oot] - tc, y[oot], 1)
        hs.append((t[w] - tc) * 24); rs.append((y[w] / np.polyval(cf, t[w] - tc) - 1) * 1e6)
    if not hs:                                  # no transit in this sector
        ax.set_title(f"Sector {sec} ({prod}): no transit covered", fontsize=9)
        continue
    h, r = np.concatenate(hs), np.concatenate(rs)
    edges = np.arange(-3 * T14 * 24, 3 * T14 * 24 + 0.01, 0.25)
    mid = (edges[1:] + edges[:-1]) / 2
    idx = np.digitize(h, edges) - 1
    med = np.array([np.mean(r[idx == i]) if np.sum(idx == i) > 2 else np.nan for i in range(len(mid))])
    err = np.array([np.std(r[idx == i]) / np.sqrt(np.sum(idx == i)) if np.sum(idx == i) > 2 else np.nan
                    for i in range(len(mid))])
    inn = np.abs(h) < 0.35 * T14 * 24
    ax.plot(h, r, ".", ms=1.5, color="0.8", zorder=0)
    ax.errorbar(mid, med, err, fmt="o", ms=3.5, color="C0", capsize=0)
    box = np.where(np.abs(mid) < T14 * 12, -DEPTH * 1e6, 0)
    ax.step(mid, box, where="mid", color="C3", lw=1.2)
    ax.set_title(f"Sector {sec} ({prod}): {len(hs)} transits, {-np.mean(r[inn]):.0f} ± "
                 f"{np.std(r[inn]) / np.sqrt(inn.sum()):.0f} ppm", fontsize=9)
    ax.axhline(0, color="k", lw=0.5)
for ax in axes.flat[len(lcs):]:                 # panels with no sector
    ax.set_visible(False)
for ax in axes[-1]:
    ax.set_xlabel("hours from mid-transit")
for ax in axes[:, 0]:
    ax.set_ylabel("flux − 1 (ppm)")
axes[0, 0].set_ylim(-3.5 * DEPTH * 1e6, 2.5 * DEPTH * 1e6)
fig.suptitle(f"TIC {TIC}: P = {P:.6f} d, folded per sector (blue: 15-min bins; "
             f"red: {DEPTH * 1e6:.0f} ppm, {T14 * 24:.2f} h)")
fig.tight_layout()
fig.savefig(OUT / "sector_folds.png", dpi=100)
print("written", OUT / "sector_folds.png")
