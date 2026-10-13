#!/usr/bin/env python3
"""
pixel_check.py - is the dip really on the target? Pixel-level checks on TESS
full-frame-image cutouts, for candidates that survive the catalogue checks.

Usage:
    python pixel_check.py results/<name>/enriched.csv               # default selection
    python pixel_check.py results/<name>/enriched.csv --target "TIC 151412755" [--candidate 1]

Default selection: rows with status 'clean', or flagged double_period_variable
with enough transits (not few_transits or single_event_dominated - with few
transits the odd/even test fires on noise), and no nearby_variable_period_match
(already explained). Any star can be checked with --target. Works on
candidates_all.csv too.

For each candidate and each TESScut sector:
  1. Aperture test: the depth in the target's brightest pixel versus a 3x3
     box. If the signal comes from a neighbour, more of it falls in the larger
     box and the depth grows.
  2. Difference image: mean out-of-dip minus mean in-dip image (per-pixel slow
     trends removed). Its centroid, converted to sky coordinates, marks the
     source of the dip; sectors are combined into one position and offset.

  3. Source star: when the dimming is off-target, the source star's own light
     curve is searched for its strongest periodicity (Lomb-Scargle), and its
     dip at the candidate's ephemeris is measured. A deep dip there means the
     source is an eclipsing binary at that period; none means the candidate
     lined up with some other variation of the source by chance.

Writes <out>/pixel_check.csv (one row per candidate and sector),
<out>/pixel_summary.csv (one verdict per candidate), and a PNG per candidate
and sector.
"""
import argparse
import logging
import re
import time
import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

warnings.filterwarnings("ignore")
import lightkurve as lk                       # noqa: E402
import matplotlib                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt               # noqa: E402
from astroquery.mast import Catalogs          # noqa: E402
from mastutil import retry, set_mast_timeout  # noqa: E402
set_mast_timeout()                    # give up on a stuck MAST query after 90 s, not 600
from astroquery.vizier import Vizier          # noqa: E402
from astropy.timeseries import LombScargle   # noqa: E402
from scipy.ndimage import median_filter       # noqa: E402

import find_planets as fp                     # noqa: E402  (light-curve download/cleaning)

# ------------------------------------------------------------------ settings
CUTOUT_PX = 11                  # cutout size (pixels; TESS pixels are 21")
PIXEL_ARCSEC = 21.0
TREND_DAYS = 3.0                # per-pixel running-median window (removes scattered light)
MAX_BIN_MIN = 30.0              # bin frames to at most this, and at most duration/6
SKY_FRAC = 0.3                  # faintest fraction of pixels used as sky in each frame
SIG_PIX = 3.0                   # pixels above this many sigma make the centroid
GROWTH_LIMIT = 1.3              # 3x3 depth / 1-pixel depth above this: off-target
OFFSET_LIMIT_PX = 0.5           # combined centroid further than this: off-target
MIN_SECTORS_FIRM = 3            # fewer centroided sectors: the verdict is only tentative
NEIGHBOUR_DMAG = 1.0            # plots: mark catalogued stars up to this much fainter than
MAX_MARKED = 8                  # the target (at most this many), and all known variables
SOURCE_TMAG_MARGIN = 3.0        # the source is the nearest star at most this much fainter
LS_PERIODS = (0.1, 15.0)        # days; periodogram range for the source's own variability
LS_EDGE_DAYS = 1.0              # drop this much from each segment edge before the periodogram
RATIOS = {1.0: "P", 2.0: "2P", 0.5: "P/2", 1.5: "3P/2", 2 / 3: "2P/3", 3.0: "3P", 1 / 3: "P/3"}
VARIABLE_CATALOGUES = (("VSX", "B/vsx/vsx", "Name", "Period"),
                       ("Gaia-EB", "I/358/veb", "Source", "Freq"),
                       ("ZTF", "J/ApJS/249/18", "ID", "Per"))
DOWNLOAD_TRIES = 4              # the cutout service throttles ("Slow Down")


# ------------------------------------------------------------------ helpers
def select_rows(df, target, candidate):
    if target:
        df = df[df["target"] == target]
        if candidate:
            df = df[df["candidate"] == candidate]
        return df
    flags = df["flags"].fillna("")
    keep = (df["status"] == "clean") | ((df["status"] == "flagged") &
                                        flags.str.contains("double_period_variable") &
                                        ~flags.str.contains("few_transits|single_event_dominated"))
    if "enrich_flags" in df:
        keep &= ~df["enrich_flags"].fillna("").str.contains("nearby_variable_period_match")
    return df[keep]


def ephemeris(r):
    """(P, T0, D, half-width of the dip) in days. Always the detected period P,
    even for double-period variables: locating the source doesn't need the two
    minima separated, and folding at 2P would leave most sectors without a dip."""
    P, T0, D = float(r["period_d"]), float(r["t0_btjd"]), float(r["duration_h"]) / 24
    W = float(r["mask_width_h"]) / 24 if "mask_width_h" in r and pd.notna(r["mask_width_h"]) else D
    return P, T0, D, max(0.5 * W, 0.5 * D)


def download(row):
    for attempt in range(1, DOWNLOAD_TRIES + 1):
        try:
            return row.download(cutout_size=CUTOUT_PX, quality_bitmask="hard")
        except Exception as e:
            logging.warning(f"    cutout attempt {attempt} failed ({type(e).__name__}); "
                            f"waiting {60 * attempt} s")
            time.sleep(60 * attempt)
    return None


def analyse(tpf, P, T0, D, half):
    """Aperture depths and difference-image centroid for one sector, or None."""
    t = tpf.time.value
    F = np.asarray(tpf.flux.value, dtype=float)
    ok = np.isfinite(t) & np.all(np.isfinite(F), axis=(1, 2))
    t, F = t[ok], F[ok]
    if len(t) < 100:
        return None
    flat = F.reshape(len(F), -1)
    F = F - np.sort(flat, axis=1)[:, : int(SKY_FRAC * flat.shape[1])].mean(axis=1)[:, None, None]
    # bin, then remove slow trends pixel by pixel
    cad = np.median(np.diff(t))
    w = max(cad, min(MAX_BIN_MIN / 1440, D / 6))     # a single dip spans several bins
    idx = np.floor((t - t.min()) / w).astype(int)
    _, inv, n = np.unique(idx, return_inverse=True, return_counts=True)
    tb = np.bincount(inv, t) / n
    ny, nx = F.shape[1:]
    Fb = np.stack([np.bincount(inv, F[:, i, j]) / n for i in range(ny) for j in range(nx)],
                  axis=1).reshape(len(tb), ny, nx)
    trend = median_filter(Fb, size=(int(TREND_DAYS / w) | 1, 1, 1), mode="nearest")
    level = np.median(Fb, axis=0)
    R = Fb - trend + level
    dt = (tb - T0 + 0.5 * P) % P - 0.5 * P
    inn = np.abs(dt) < 0.4 * D                                   # dip core
    out = (np.abs(dt) > half + 0.5 * D) & (np.abs(dt) < half + 2.5 * D)   # local baseline
    if inn.sum() < 3 or out.sum() < 6:
        return "no dip in this sector" if inn.sum() == 0 else "too little data around the dip"
    mo, mi = R[out].mean(axis=0), R[inn].mean(axis=0)
    diff = mo - mi                                               # >0 where it dims
    err = R[out].std(axis=0) * np.sqrt(1 / inn.sum() + 1 / out.sum())

    wcs = tpf.wcs
    tx, ty = wcs.world_to_pixel_values(tpf.ra, tpf.dec)
    ix, iy = int(round(float(tx))), int(round(float(ty)))
    ix, iy = min(max(ix, 1), nx - 2), min(max(iy, 1), ny - 2)
    box = (slice(iy - 1, iy + 2), slice(ix - 1, ix + 2))
    d1, e1 = diff[iy, ix] / mo[iy, ix], err[iy, ix] / mo[iy, ix]
    d3 = diff[box].sum() / mo[box].sum()
    e3 = np.sqrt((err[box] ** 2).sum()) / mo[box].sum()

    z = diff / err
    sel = z > SIG_PIX
    res = dict(n_in=int(inn.sum()), depth_1px=d1, err_1px=e1, depth_3x3=d3, err_3x3=e3,
               peak_sigma=float(z.max()), ra=np.nan, dec=np.nan,
               img=level, diff=diff, tx=float(tx), ty=float(ty), cx=np.nan, cy=np.nan, wcs=wcs)
    if sel.any():
        yy, xx = np.indices(diff.shape)
        cx = (diff[sel] * xx[sel]).sum() / diff[sel].sum()
        cy = (diff[sel] * yy[sel]).sum() / diff[sel].sum()
        c = wcs.pixel_to_world(cx, cy)
        res.update(cx=cx, cy=cy, ra=c.ra.deg, dec=c.dec.deg)
    return res


def plot(path, title, r, field):
    """Median and difference images, with the target, the dimming centroid,
    brighter catalogued neighbours (white circles, TESS mag) and known
    variables (orange diamonds) marked."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    ny, nx = r["img"].shape
    marks = []
    for kind, items in (("star", field["stars"]), ("var", field["variables"])):
        for it in items:
            x, y = r["wcs"].world_to_pixel_values(it["ra"], it["dec"])
            if -0.5 <= x <= nx - 0.5 and -0.5 <= y <= ny - 0.5:
                marks.append((kind, float(x), float(y), it["label"]))
    for ax, img, label in ((axes[0], r["img"], "Median image"),
                           (axes[1], r["diff"], "Out-of-dip minus in-dip (dimming)")):
        im = ax.imshow(img, origin="lower", cmap="viridis")
        fig.colorbar(im, ax=ax, shrink=0.8)
        shown = set()
        for kind, x, y, text in marks:
            style = (dict(marker="o", mfc="none", mec="w", ms=11, mew=1.5, ls="none")
                     if kind == "star" else
                     dict(marker="D", mfc="none", mec="orange", ms=13, mew=2, ls="none"))
            name = "brighter stars (T mag)" if kind == "star" else "known variables"
            ax.plot(x, y, label=None if name in shown else name, **style)
            shown.add(name)
            ax.annotate(text, (x, y), color="w" if kind == "star" else "orange", fontsize=7,
                        xytext=(5, 5 if kind == "star" else -9), textcoords="offset points")
        ax.plot(r["tx"], r["ty"], "rx", ms=12, mew=2, label="target")
        if np.isfinite(r["cx"]):
            ax.plot(r["cx"], r["cy"], "m+", ms=14, mew=2, label="dimming centroid")
        ax.set_title(label)
    axes[0].legend(loc="lower left", fontsize=7, framealpha=0.6)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def field_objects(coord):
    """Catalogued objects across the cutout, for marking on the plots: the
    target's TESS magnitude, brighter neighbours and known variables."""
    radius = CUTOUT_PX / 2 * np.sqrt(2) * PIXEL_ARCSEC * u.arcsec     # cutout corners
    field = {"stars": [], "variables": []}
    try:
        tab = Catalogs.query_region(coord, radius=radius, catalog="TIC").to_pandas()
        tab = tab[np.isfinite(tab["Tmag"])].sort_values("dstArcSec")
        t_target = tab.iloc[0]["Tmag"]
        field["tmag"] = float(t_target)
        near = tab[(tab["dstArcSec"] > 2) & (tab["Tmag"] <= t_target + NEIGHBOUR_DMAG)]
        field["stars"] = [dict(ra=s["ra"], dec=s["dec"], label=f"{s['Tmag']:.1f}")
                          for _, s in near.sort_values("Tmag").head(MAX_MARKED).iterrows()]
    except Exception as e:
        logging.warning(f"  TIC field query failed ({type(e).__name__}); neighbours not marked")
    v = Vizier(columns=["**"], row_limit=-1)
    seen = []
    for label, cat, name, per in VARIABLE_CATALOGUES:
        try:
            res = v.query_region(coord, radius=radius, catalog=cat)
        except Exception:
            continue
        tabs = [x for x in res.values() if name in x.colnames and per in x.colnames]
        if not tabs:
            continue
        tab = tabs[0]
        ra_col = next(c for c in ("RAJ2000", "RA_ICRS") if c in tab.colnames)
        dec_col = ra_col.replace("RA", "DE")
        for r in tab:
            ra, dec = float(r[ra_col]), float(r[dec_col])
            if any(abs(ra - a) * np.cos(np.radians(dec)) * 3600 < 3 and abs(dec - b) * 3600 < 3
                   for a, b in seen):               # same star in another catalogue
                continue
            seen.append((ra, dec))
            p = float(r[per])
            p = 1 / p if per == "Freq" and p > 0 else p
            kind = str(r["Type"]) if "Type" in r.colnames else "EB"
            field["variables"].append(dict(ra=ra, dec=dec, label=kind + (f" {p:.4g}d" if p > 0 else "")))
    return field


def source_star(coord, offset, tmag=np.nan):
    """What sits at the dimming position: the nearest TIC stars and any known
    variables (VSX, Gaia DR3 eclipsing binaries, ZTF), each with its distance
    from that position. The search radius stays well short of the target
    (offset = its distance), so the target itself is never named."""
    radius = min(30.0, 0.6 * offset) * u.arcsec
    parts, source = [], None
    try:
        tab = Catalogs.query_region(coord, radius=radius, catalog="TIC").to_pandas()
        tab = tab[np.isfinite(tab["Tmag"])].sort_values("dstArcSec")
        bright = tab if not np.isfinite(tmag) else tab[tab["Tmag"] <= tmag + SOURCE_TMAG_MARGIN]
        if len(bright) or len(tab):
            source = str((bright if len(bright) else tab).iloc[0]["ID"])
        parts += [f"TIC {s['ID']} (T={s['Tmag']:.1f}, {s['dstArcSec']:.0f}\")"
                  for _, s in tab.head(3).iterrows()]
    except Exception:
        pass
    v = Vizier(columns=["**", "+_r"], row_limit=-1)
    found = []
    for label, cat, name, per in VARIABLE_CATALOGUES:
        try:
            res = v.query_region(coord, radius=radius, catalog=cat)
        except Exception:
            continue
        tabs = [x for x in res.values() if name in x.colnames and per in x.colnames]
        for r in (tabs[0] if tabs else []):
            ra_col = next(c for c in ("RAJ2000", "RA_ICRS") if c in r.colnames)
            sep = SkyCoord(float(r[ra_col]) * u.deg, float(r[ra_col.replace("RA", "DE")]) * u.deg
                           ).separation(coord).arcsec
            p = float(r[per])
            p = 1 / p if per == "Freq" and p > 0 else p
            kind = str(r["Type"]) if "Type" in r.colnames else "EB"
            found.append(f"{label} {r[name]} ({kind}" + (f", P={p:.5f} d" if p > 0 else "")
                         + f", {sep:.0f}\")")
    if found:
        parts.append("known variable(s): " + "; ".join(found))
    return "; ".join(parts), source


def remove_periodic(T, Y, Pv, P, T0, D, harmonics=4):
    """Subtract a Fourier series at period Pv (plus a quadratic) from each data
    segment, fitted while ignoring the candidate's dip times. Returns (t, y)."""
    out_t, out_y = [], []
    for ts, ys in zip(T, Y):
        ph = 2 * np.pi * ts / Pv
        x = ts - ts.mean()
        A = np.column_stack([np.ones_like(ts), x, x ** 2] +
                            [f(k * ph) for k in range(1, harmonics + 1) for f in (np.sin, np.cos)])
        dt = (ts - T0 + 0.5 * P) % P - 0.5 * P
        fit = np.abs(dt) > D                              # leave the dips out of the fit
        if fit.sum() < A.shape[1] * 3:
            continue
        coef, *_ = np.linalg.lstsq(A[fit], ys[fit], rcond=None)
        out_t.append(ts)
        out_y.append(ys - A @ coef)
    if not out_t:
        return np.array([]), np.array([])
    return np.concatenate(out_t), np.concatenate(out_y)


def ephemeris_dip(t, y, P, T0, D, n_random=60):
    """Mean flux just outside minus inside the candidate's dips, and its error
    from the same measurement at n_random other times (so slow, correlated
    variability is included). NaNs if there is no data in the dips."""
    def dip(t0):
        dt = (t - t0 + 0.5 * P) % P - 0.5 * P
        inn = np.abs(dt) < 0.5 * D
        near = (np.abs(dt) > 0.75 * D) & (np.abs(dt) < 2.5 * D)
        if inn.sum() < 3 or near.sum() < 6:
            return np.nan
        return np.mean(y[near]) - np.mean(y[inn])
    d = dip(T0)
    if not np.isfinite(d):
        return np.nan, np.nan
    rng = np.random.default_rng(0)
    shifts = rng.uniform(3 * D, P - 3 * D, n_random) if P > 6 * D else rng.uniform(0, P, n_random)
    ctrl = np.array([dip(T0 + s) for s in shifts])
    ctrl = ctrl[np.isfinite(ctrl)]
    err = np.std(ctrl) if len(ctrl) >= 10 else np.nan
    return d, err


def source_variability(tic, P, T0, D):
    """The source star's own light curve: its strongest periodicity and its dip
    at the candidate's ephemeris. Text, or '' if no light curve is available."""
    name = f"TIC {tic}"
    sr = retry(lambda: lk.search_lightcurve(name, author=fp.AUTHORS), f"{name}: MAST search")
    if len(sr) == 0:
        return f"source {name}: no light curve available"
    rows, _ = fp.choose_products(sr)
    lcc = fp.download(sr[rows]) if rows else None
    lcs = [c for c in (fp.clean(lc) for lc in (lcc or [])) if c is not None]
    if not lcs:
        return f"source {name}: no usable light curve"
    # periodogram: per segment, drop the edges (ramps) and normalise
    T, Y = [], []
    for c in lcs:
        t, y = fp.arr(c.time), fp.arr(c.flux)
        for seg in np.split(np.arange(len(t)), np.flatnonzero(np.diff(t) > 1.0) + 1):
            ts, ys = t[seg], y[seg]
            ok = (ts - ts[0] > LS_EDGE_DAYS) & (ts[-1] - ts > LS_EDGE_DAYS)
            if ok.sum() > 50:
                T.append(ts[ok])
                Y.append(ys[ok] / np.median(ys[ok]) - 1)
    text = f"source {name}:"
    if T:
        t, y = np.concatenate(T), np.concatenate(Y)
        freq, power = LombScargle(t, y).autopower(minimum_frequency=1 / LS_PERIODS[1],
                                                 maximum_frequency=1 / LS_PERIODS[0],
                                                 samples_per_peak=10)
        Pv = 1 / freq[np.argmax(power)]
        ph = (t / Pv) % 1
        prof = [np.median(y[(ph >= a) & (ph < a + 0.05)]) for a in np.arange(0, 1, 0.05)
                if np.any((ph >= a) & (ph < a + 0.05))]
        amp = (max(prof) - min(prof)) / 2
        rel = next((f"= {lab} of the candidate's period" for f, lab in RATIOS.items()
                    if abs(Pv / (f * P) - 1) < 0.01), f"unrelated to the candidate's period "
                                                       f"(ratio {P / Pv:.2f})")
        text += f" strongest periodicity {Pv:.4f} d, amplitude +-{amp * 100:.1f}% ({rel});"
    # dip at the candidate's ephemeris
    related = T and not rel.startswith("unrelated")
    if T and not related:
        # The source's own variability is something else: remove it first, or
        # its minima (which the candidate may have lined up with) read as a dip.
        t, y = remove_periodic(T, Y, Pv, P, T0, D)
    else:
        parts = [fp.flatten_sector(c, [(P, T0, D)]) for c in lcs]
        t, y = (np.concatenate([p[i] for p in parts]) for i in (0, 1))
        y = y - 1
    d, err = ephemeris_dip(t, y, P, T0, D)
    if np.isfinite(d):
        text += (f" at the candidate's dip times it dims by {d * 1e6:.0f} +- {err * 1e6:.0f} ppm "
                 f"({d / err:.1f} sigma vs random times"
                 + (f", after removing its {Pv:.4f} d variation)" if T and not related else ")"))
    else:
        text += " no data at the candidate's dip times"
    return text


def check(r, outdir):
    name, k = r["target"], int(r["candidate"])
    P, T0, D, half = ephemeris(r)
    logging.info(f"{name} #{k}: P = {P:.6f} d, T0 = {T0:.4f}, dur = {D * 24:.2f} h")
    sr = retry(lambda: lk.search_tesscut(name), f"{name}: TESScut search")
    safe = re.sub(r"[^\w.-]+", "_", name)
    rows, tgt, field = [], None, None
    for row in sr:
        sector = int(row.table["sequence_number"][0])
        tpf = download(row)
        if tpf is None:
            logging.warning(f"  sector {sector}: download failed, skipped")
            continue
        tgt = SkyCoord(tpf.ra * u.deg, tpf.dec * u.deg)
        if field is None:
            field = field_objects(tgt)
        a = analyse(tpf, P, T0, D, half)
        if a is None or isinstance(a, str):
            logging.info(f"  sector {sector}: {a or 'too little good data'}, skipped")
            continue
        off = (SkyCoord(a["ra"] * u.deg, a["dec"] * u.deg).separation(tgt).arcsec
               if np.isfinite(a["ra"]) else np.nan)
        logging.info(f"  sector {sector}: depth 1px {a['depth_1px'] * 1e6:.0f}+-{a['err_1px'] * 1e6:.0f}"
                     f" ppm, 3x3 {a['depth_3x3'] * 1e6:.0f}+-{a['err_3x3'] * 1e6:.0f} ppm; "
                     f"difference peak {a['peak_sigma']:.0f} sigma, centroid {off:.0f}\" from target")
        plot(outdir / f"{safe}_c{k}_s{sector}.png",
             f"{name} #{k} sector {sector}: P = {P:.5f} d", a, field)
        rows.append(dict(target=name, candidate=k, sector=sector, period_d=P, n_in=a["n_in"],
                         depth_1px_ppm=a["depth_1px"] * 1e6, err_1px_ppm=a["err_1px"] * 1e6,
                         depth_3x3_ppm=a["depth_3x3"] * 1e6, err_3x3_ppm=a["err_3x3"] * 1e6,
                         peak_sigma=a["peak_sigma"], centroid_ra=a["ra"], centroid_dec=a["dec"],
                         offset_arcsec=off))
    if not rows:
        logging.info("  -> no usable sectors")
        return rows, dict(target=name, candidate=k, verdict="no usable sectors", source="")
    df = pd.DataFrame(rows)
    # aperture test: weighted mean ratio over sectors with a detected 1-pixel dip
    good = df[df.depth_1px_ppm > 3 * df.err_1px_ppm]
    verdict, source_text = [], ""
    if len(good):
        wts = (good.depth_1px_ppm / good.err_1px_ppm) ** 2
        growth = np.sum(wts * good.depth_3x3_ppm / good.depth_1px_ppm) / wts.sum()
        verdict.append(f"3x3/1px depth = {growth:.2f}"
                       + (" (grows: off-target?)" if growth > GROWTH_LIMIT else ""))
    else:
        verdict.append("no sector detects the dip in the target's own pixel")
    # combined centroid, weighted by significance
    cen = df[np.isfinite(df.centroid_ra)]
    if len(cen):
        wts = cen.peak_sigma ** 2
        ra = np.sum(wts * cen.centroid_ra) / wts.sum()
        dec = np.sum(wts * cen.centroid_dec) / wts.sum()
        c = SkyCoord(ra * u.deg, dec * u.deg)
        seps = np.array([SkyCoord(a * u.deg, b * u.deg).separation(c).arcsec
                         for a, b in zip(cen.centroid_ra, cen.centroid_dec)])
        scatter = np.sqrt(np.mean(seps ** 2))          # rms distance from the mean position
        off = c.separation(tgt).arcsec
        verdict.append(f"dimming centre {off:.0f}\" from target (sector scatter {scatter:.0f}\", "
                       f"{len(cen)} sectors)")
        firm = len(cen) >= MIN_SECTORS_FIRM
        if off > OFFSET_LIMIT_PX * PIXEL_ARCSEC and off > 2 * scatter / np.sqrt(len(cen)):
            star, tic = source_star(c, off, (field or {}).get("tmag", np.nan))
            verdict.append(("OFF-TARGET" if firm else "possibly off-target (too few sectors)")
                           + f" at RA {ra:.4f} Dec {dec:.4f}" + (f": {star}" if star else ""))
            if tic:
                try:
                    source_text = source_variability(tic, P, T0, D)
                except Exception as e:
                    source_text = f"source TIC {tic}: light-curve check failed ({type(e).__name__})"
        else:
            verdict.append("consistent with the target"
                           + ("" if firm else " (too few sectors to be sure)"))
    logging.info("  -> " + "; ".join(verdict))
    if source_text:
        logging.info("     " + source_text)
    return rows, dict(target=name, candidate=k, verdict="; ".join(verdict), source=source_text)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", help="enriched.csv or candidates_all.csv")
    ap.add_argument("--target", help="check this star only")
    ap.add_argument("--candidate", type=int, help="with --target: this candidate number only")
    ap.add_argument("--out", help="output folder (default: pixel_check beside the input)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("lightkurve").setLevel(logging.WARNING)   # quality-mask chatter
    src = Path(args.candidates)
    outdir = Path(args.out) if args.out else src.parent / "pixel_check"
    outdir.mkdir(parents=True, exist_ok=True)
    df = select_rows(pd.read_csv(src), args.target, args.candidate)
    logging.info(f"{len(df)} candidate(s) to check")
    rows, summary = [], []
    for _, r in df.iterrows():
        try:
            sector_rows, verdict = check(r, outdir)
            rows += sector_rows
            summary.append(verdict)
        except Exception:
            logging.exception(f"  {r['target']} #{r['candidate']} failed")
    if rows:
        pd.DataFrame(rows).to_csv(outdir / "pixel_check.csv", index=False)
    if summary:
        pd.DataFrame(summary).to_csv(outdir / "pixel_summary.csv", index=False)
        logging.info(f"wrote {len(rows)} sector rows to {outdir / 'pixel_check.csv'}")


if __name__ == "__main__":
    main()
