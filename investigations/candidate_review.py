"""Close look at one candidate: SPOC DV check, per-event fits and timing, and a
centroid test against the brightest neighbours.

Usage: python investigations/candidate_review.py "TIC 298514156" P T0_BTJD DURATION_H OUTDIR
"""
import sys
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path as _Path
ROOT = str(_Path(__file__).resolve().parents[1])   # the repository
sys.path.insert(0, ROOT)
from pathlib import Path

import astropy.units as u
import lightkurve as lk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.mast import Catalogs, Observations
from scipy.optimize import least_squares

import find_planets as fp

name, P0, T00, D0 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) / 24
out = Path(sys.argv[5]); out.mkdir(parents=True, exist_ok=True)
tic = name.split()[-1]

# ---- 1. did SPOC's pipeline flag anything (DV products)?
obs = Observations.query_criteria(obs_collection="TESS", target_name=tic)
prods = Observations.get_product_list(obs) if len(obs) else []
dv = [p for p in prods if str(p["productSubGroupDescription"]).upper().startswith("DV")]
print(f"1. SPOC data-validation products on MAST: {len(dv)}"
      + ("" if not dv else " -> " + ", ".join(sorted({str(p['productFilename']) for p in dv})[:6])))

# ---- 2. events: local quadratic baseline, common trapezoid, individual times
sr = lk.search_lightcurve(name, author=fp.AUTHORS)
rows, label = fp.choose_products(sr)
lcs = [c for c in (fp.clean(lc) for lc in sr[rows].download_all()) if c is not None]
WIN, CORE = max(0.75, 4 * D0), 1.5 * D0
events = []
for c in lcs:
    t, y, e = fp.arr(c.time), fp.arr(c.flux), fp.arr(c.flux_err)
    for k in np.arange(np.ceil((t.min() - 1 - T00) / P0), np.floor((t.max() + 1 - T00) / P0) + 1):
        tc = T00 + k * P0
        sel = np.abs(t - tc) < WIN
        oot = sel & (np.abs(t - tc) > CORE)
        if np.sum(np.abs(t - tc) < D0 / 2) < 5 or np.sum(oot & (t < tc)) < 20 or np.sum(oot & (t > tc)) < 20:
            continue
        cf = np.polyfit(t[oot] - tc, y[oot], 2)
        base = np.polyval(cf, t[sel] - tc)
        events.append(dict(sector=c.meta["SECTOR"], epoch=int(k), t=t[sel], y=y[sel] / base,
                           e=e[sel] / base, cad=np.median(np.diff(t))))
print(f"2. data: {label}; {len(events)} events with data on both sides")


def trap(dt, depth, t14, tau):
    a = np.abs(dt); half = t14 / 2
    return np.where(a <= half - tau, depth, np.where((a > half - tau) & (a < half),
                    depth * (half - a) / max(tau, 1e-6), 0.0))


def model(ev, tc, depth, t14, tau):
    sub = (np.arange(5) - 2) / 5 * ev["cad"]
    return 1 - np.mean([trap(ev["t"] + s - tc, depth, t14, tau) for s in sub], axis=0)


def resid(p):
    T0, P, depth, t14, tau = p
    return np.concatenate([(ev["y"] - model(ev, T0 + ev["epoch"] * P, depth, t14, tau)) / ev["e"] for ev in events])


fit = least_squares(resid, [T00, P0, 3e-4, D0, 0.25 * D0],
                    bounds=([T00 - D0, P0 - 0.01, 0, 0.3 * D0, 0.02 * D0], [T00 + D0, P0 + 0.01, 0.01, 3 * D0, 1.5 * D0]))
chi2 = np.sum(fit.fun ** 2); dof = len(fit.fun) - 5
err = np.sqrt(np.diag(np.linalg.inv(fit.jac.T @ fit.jac) * chi2 / dof))
T0, P, depth, t14, tau = fit.x
print(f"   global fit (chi2/dof {chi2 / dof:.2f}): P = {P:.5f} +- {err[1]:.5f} d, T0 = BTJD {T0:.4f} +- {err[0]:.4f}, "
      f"depth {depth * 1e6:.0f} +- {err[2] * 1e6:.0f} ppm, T14 {t14 * 24:.2f} +- {err[3] * 24:.2f} h, ingress {tau * 24:.2f} h")
fig, axes = plt.subplots(len(events), 1, figsize=(9, 1.9 * len(events)), sharex=True)
for ax, ev in zip(np.atleast_1d(axes), events):
    tcp = T0 + ev["epoch"] * P
    grid = tcp + np.linspace(-1.5, 1.5, 361) * t14
    c2 = np.array([np.sum(((ev["y"] - model(ev, g, depth, t14, tau)) / ev["e"]) ** 2) for g in grid])
    flat = np.sum(((ev["y"] - 1) / ev["e"]) ** 2)
    i = int(np.argmin(c2)); red = c2[i] / max(1, len(ev["y"]) - 1)
    inside = grid[c2 <= c2[i] + red]
    sig = (inside.max() - inside.min()) / 2
    at_pred = np.sum(((ev["y"] - model(ev, tcp, depth, t14, tau)) / ev["e"]) ** 2)
    print(f"   S{ev['sector']} E{ev['epoch']}: O-C {(grid[i] - tcp) * 1440:+5.0f} +- {sig * 1440:3.0f} min; "
          f"transit vs none at predicted time: delta-chi2 {flat - at_pred:6.1f} (rescaled {(flat - at_pred) / red:6.1f})")
    dt = (ev["t"] - tcp) * 24
    b = np.arange(dt.min(), dt.max(), 0.33)
    bm = [np.median(ev["y"][(dt >= a) & (dt < a + 0.33)]) for a in b]
    ax.plot(dt, (ev["y"] - 1) * 1e6, ".", ms=1.2, color="0.75")
    ax.plot(b + 0.165, (np.array(bm) - 1) * 1e6, "o-", ms=2.5, color="C0", lw=1)
    ax.plot(dt, (model(ev, tcp, depth, t14, tau) - 1) * 1e6, color="C3", lw=1.5)
    ax.set_ylim(-1500, 1000); ax.text(0.01, 0.85, f"S{ev['sector']} E{ev['epoch']}", transform=ax.transAxes, fontsize=8)
np.atleast_1d(axes)[-1].set_xlabel("hours from predicted mid-transit")
fig.suptitle(f"{name}: events (grey 2-min / FFI, blue 20-min medians, red common fit)", fontsize=10)
fig.tight_layout(); fig.savefig(out / "events.png", dpi=100)

# ---- 3. centroid test (SPOC 2-min), with predictions for each bright neighbour
star = Catalogs.query_criteria(catalog="TIC", ID=int(tic))[0]
tgt = SkyCoord(float(star["ra"]) * u.deg, float(star["dec"]) * u.deg)
nb = Catalogs.query_region(tgt, radius=2 * u.arcmin, catalog="TIC").to_pandas()
nb = nb[(nb["ID"].astype(str) != tic) & (nb["Tmag"] < float(star["Tmag"]) + 5)].sort_values("Tmag")
spoc = lk.search_lightcurve(name, author="SPOC", exptime=120)
dcol, drow, wts, pred = [], [], [], {}
for r in spoc:
    raw = r.download(quality_bitmask="default")
    tp = lk.search_targetpixelfile(name, author="SPOC", exptime=120, sector=raw.meta["SECTOR"]).download(quality_bitmask="default")
    t = np.asarray(raw.time.value); cx = np.asarray(raw.mom_centr1.value, float); cy = np.asarray(raw.mom_centr2.value, float)
    ok = np.isfinite(t) & np.isfinite(cx) & np.isfinite(cy)
    for k in np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1):
        tc = T0 + k * P
        inn = ok & (np.abs(t - tc) < 0.4 * t14)
        oot = ok & (np.abs(t - tc) > 0.75 * t14) & (np.abs(t - tc) < max(0.5, 4 * t14))
        if inn.sum() < 10 or oot.sum() < 40:
            continue
        for v, store in ((cx, dcol), (cy, drow)):
            cf = np.polyfit(t[oot] - tc, v[oot], 1)
            dv = v - np.polyval(cf, t - tc)
            store.append((np.mean(dv[inn]), np.std(dv[oot]) / np.sqrt(inn.sum())))
        # photocentre and neighbour positions (pixels) in this sector
        xc, yc = np.median(cx[oot]), np.median(cy[oot])
        for _, s in nb.head(6).iterrows():
            x, y = tp.wcs.world_to_pixel_values(s["ra"], s["dec"])
            x, y = x + tp.column, y + tp.row
            pred.setdefault(s["ID"], []).append((-depth * (x - xc), -depth * (y - yc), s["Tmag"], s["dstArcSec"]))
        x, y = tp.wcs.world_to_pixel_values(tgt.ra.deg, tgt.dec.deg)
        pred.setdefault("target", []).append((-depth * (x + tp.column - xc), -depth * (y + tp.row - yc),
                                              float(star["Tmag"]), 0.0))
        print(f"   centroid S{raw.meta['SECTOR']} E{int(k)}: col {dcol[-1][0] * 21000:+6.1f} +- {dcol[-1][1] * 21000:4.1f}, "
              f"row {drow[-1][0] * 21000:+6.1f} +- {drow[-1][1] * 21000:4.1f} mas")
if dcol:
    def wmean(v):
        a = np.array(v); w = 1 / a[:, 1] ** 2
        return np.sum(a[:, 0] * w) / w.sum(), 1 / np.sqrt(w.sum())
    (mx, ex), (my, ey) = wmean(dcol), wmean(drow)
    print(f"3. measured in-transit centroid shift: col {mx * 21000:+.1f} +- {ex * 21000:.1f}, "
          f"row {my * 21000:+.1f} +- {ey * 21000:.1f} mas")
    print("   predicted shift if the dip came from (same depth, same events; mas):")
    for key, v in pred.items():
        a = np.array([x[:2] for x in v]).mean(axis=0) * 21000
        tm, sep = v[0][2], v[0][3]
        chi = np.hypot((mx * 21000 - a[0]) / (ex * 21000), (my * 21000 - a[1]) / (ey * 21000))
        lab = "the target" if key == "target" else f"TIC {key} (T={tm:.1f}, {sep:.0f}\")"
        print(f"     {lab:34s} col {a[0]:+6.1f}, row {a[1]:+6.1f}  -> measured differs by {chi:4.1f} sigma")
