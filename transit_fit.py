#!/usr/bin/env python3
"""
transit_fit.py - physical transit fit and stellar-density check for clean candidates.

Usage:
    python transit_fit.py results/<name>/enriched.csv            # clean candidates
    python transit_fit.py results/<name>/enriched.csv --target "TIC 302453503" [--candidate 1]
    python transit_fit.py results/<name>/candidates_all.csv --all  # every candidate (calibration)
    options: --out DIR (default: transit_fit beside the input), --steps N (most MCMC steps, default 30000;
             the chain stops earlier once converged: > 50 autocorrelation times, stable to 1%),
             --jobs N (fits run in parallel, default 2); MAST queries are retried

A transit's shape and duration fix a/R*, and with the period the star's density:
rho* = 3 pi (a/R*)^3 / (G P^2). If the dip is on the catalogued star, that density
matches the TIC's (from its mass and radius). A blend with a different star, or an
eclipsing binary, usually doesn't: the check that sank the candidate in a
Planet Hunters TESS thread (8 vs 2 g/cm3) and that TIC 302453503 passes.

For each candidate this script loads the search's light curves, normalises each
transit by a quadratic fitted outside +-1.5 durations (up to MAX_EVENTS transits,
spread across the data), and runs vetting/transit_mcmc.py (batman + emcee) in the
vetting environment (.venv-vet; see requirements-vet.txt). Limb-darkening priors
come from the star's Teff (approximate TESS-band values).

Verdict, from log10 density (fit posterior combined with the TIC uncertainty):
  consistent        the TIC density is within DENSITY_SIGMA
  density_mismatch  it is further off (flag: density_mismatch)
  grazing_density_uncertain
                    a density_mismatch on a v_shaped signal: the fit can't pin b and Rp/R*
                    for a grazing transit (TOI-7689.01 converges on b ~ 0.57 with ~14x the
                    star's density), so the density test isn't trusted, as PHT treats
                    V-shaped candidates (Eisner et al. 2021). Goes on to the FPP
  duration_too_long the fitted T14 is over DUR_MARGIN x the longest transit the TIC star
                    allows (central, circular, with the fitted Rp/R*). Takes precedence:
                    such a fit can still pass the density test through wide errors (a
                    b023 blend: 7.7 h vs 4.4 h; a b017 fit that wandered to 9.9 h).
                    Real candidates so far reach 1.16x
Only consistent fits go on to fpp_check.py.
  no_stellar_params the TIC has no mass or radius
  fit_failed
Writes <out>/transit_fit.csv; the events are kept in <out>/ for fpp_check.py.
"""
import argparse
import json
import logging
import re
import subprocess
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

import logsetup

warnings.filterwarnings("ignore")
from mastutil import set_mast_timeout  # noqa: E402
set_mast_timeout()                    # give up on a stuck MAST query after 90 s, not 600
HERE = Path(__file__).resolve().parent
VET_PY = HERE / ".venv-vet" / "bin" / "python"
MAX_EVENTS = 60             # transits fitted (spread evenly); more add little but cost time
DENSITY_SIGMA = 3.0         # density_mismatch beyond this many sigma (log10 density)
TTV_P = 0.001               # note a density_mismatch whose sector timing scatters this significantly
TRIES = 3                   # MAST queries: attempts, 60 s apart
RHO_SUN = 1.41              # g/cm3
DUR_MARGIN = 1.3            # duration_too_long above this multiple of the star's longest transit


def retry(func, what):
    """Run a MAST/network call, retrying transient failures (timeouts)."""
    import time
    for attempt in range(1, TRIES + 1):
        try:
            return func()
        except Exception as e:
            if attempt == TRIES:
                raise
            logging.warning(f"  {what} failed ({type(e).__name__}); retrying in 60 s")
            time.sleep(60)


def ld_prior(teff):
    """Approximate quadratic limb darkening in the TESS band (after Claret 2017),
    used as prior means (width 0.1)."""
    if not np.isfinite(teff):
        return 0.30, 0.30
    if teff < 3800:
        return 0.20, 0.35
    if teff < 5000:
        return 0.38, 0.22
    if teff < 6000:
        return 0.35, 0.25
    return 0.28, 0.28


def stellar_params(name):
    """TIC radius, mass (with errors) and Teff, for "TIC <n>" or a resolvable name."""
    from astroquery.mast import Catalogs
    m = re.fullmatch(r"TIC\s*(\d+)", name)
    if m:
        r = Catalogs.query_criteria(catalog="TIC", ID=int(m.group(1)))
    else:
        r = Catalogs.query_object(name, catalog="TIC", radius=0.002)
        if len(r):
            r.sort("dstArcSec")
    if not len(r):
        return None
    r = r[0]
    num = lambda k: float(r[k]) if r[k] is not None and np.isfinite(float(r[k] or np.nan)) else np.nan
    return dict(rad=num("rad"), e_rad=num("e_rad"), mass=num("mass"), e_mass=num("e_mass"), teff=num("Teff"))


def local_events(lcs, P, T0, D):
    """Each transit with data on both sides, normalised by a quadratic fitted
    outside +-1.5 D within +-max(0.75 d, 4 D)."""
    import find_planets as fp
    win, core = max(0.75, 4 * D), 1.5 * D
    ev = []
    for c in lcs:
        t, y, e = fp.arr(c.time), fp.arr(c.flux), fp.arr(c.flux_err)
        ok = np.isfinite(t) & np.isfinite(y) & np.isfinite(e)
        t, y, e = t[ok], y[ok], e[ok]
        cad = np.median(np.diff(t))
        for k in np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1):
            tc = T0 + k * P
            sel = np.abs(t - tc) < win
            oot = sel & (np.abs(t - tc) > core)
            if (np.sum(np.abs(t - tc) < D / 2) < max(3, 0.3 * D / cad) or np.sum(oot & (t < tc)) < 10
                    or np.sum(oot & (t > tc)) < 10):
                continue
            base = np.polyval(np.polyfit(t[oot] - tc, y[oot], 2), t[sel] - tc)
            near = np.abs(t[sel] - tc) < max(3 * D, 0.15)       # keep the fit's window only
            ev.append((t[sel][near], (y[sel] / base)[near], (e[sel] / base)[near],
                       np.full(near.sum(), cad), c.meta.get("SECTOR")))
    if len(ev) > MAX_EVENTS:
        ev = [ev[i] for i in np.linspace(0, len(ev) - 1, MAX_EVENTS).round().astype(int)]
    return ev


def sector_timing(ev, P, T0, D, depth):
    """Diagnostic only. Transit-timing variations (e.g. planets near resonance,
    like TOI-270 c, d) shift whole groups of transits by minutes; stacked on a
    fixed ephemeris the dip smears into a V, which the fit reads as a grazing
    transit of a low-density star. Each sector's timing offset is measured with
    a fixed trapezoid shape and tested for scatter about a straight line (a
    linear drift is just the ephemeris). Aligning sectors before the fit was
    tried and rejected: it fixed TOI-270 c, d, but for weak signals the offsets
    are mostly noise, and aligning them sharpens noise into a dip (TOI-700 e
    "corrected" by 116 min; pi Men c pushed to a mismatch). So the result is only
    reported, as a caution on a density_mismatch. Returns (events, info)."""
    from scipy.optimize import least_squares
    from scipy.stats import chi2 as chi2d
    import find_planets as fp
    t = np.concatenate([x[0] for x in ev]); y = np.concatenate([x[1] for x in ev])
    e = np.concatenate([x[2] for x in ev])
    ph = lambda tt, sh=0.0: (tt - sh - T0 + 0.5 * P) % P - 0.5 * P
    res = least_squares(lambda q: (y - (1 - fp.trapezoid_drop(ph(t), q[0], q[1], q[2] * q[1] / 2))) / e,
                        [depth, D, 0.3], bounds=([0, 0.3 * D, 0.02], [0.05, 3 * D, 1.0]))
    dep, t14, f = res.x
    shifts = np.linspace(-0.5 * D, 0.5 * D, 121)
    secs, rows = sorted({x[4] for x in ev}, key=lambda v: (v is None, v)), []
    for sct in secs:
        g = [x for x in ev if x[4] == sct]
        tt = np.concatenate([x[0] for x in g]); yy = np.concatenate([x[1] for x in g])
        ee = np.concatenate([x[2] for x in g])
        c2 = np.array([np.sum(((yy - (1 - fp.trapezoid_drop(ph(tt, sh), dep, t14, f * t14 / 2))) / ee) ** 2)
                       for sh in shifts])
        i = int(np.argmin(c2))
        red = max(1.0, c2[i] / max(1, len(tt) - 1))
        inside = shifts[c2 <= c2[i] + red]
        err = max((inside.max() - inside.min()) / 2, shifts[1] - shifts[0])
        if 0 < i < len(shifts) - 1:                       # an offset at the edge of the scan isn't measured
            rows.append((sct, float(np.mean(tt)), float(shifts[i]), float(err)))
    info = dict(ttv_p=np.nan, ttv_max_min=np.nan, ttv_sectors=len(rows))
    if len(rows) < 3:
        return ev, info
    tm, off, er = (np.array([r_[k] for r_ in rows]) for k in (1, 2, 3))
    w = 1 / er ** 2
    A = np.vstack([np.ones_like(tm), tm - tm.mean()]).T
    coef = np.linalg.lstsq(A * np.sqrt(w)[:, None], off * np.sqrt(w), rcond=None)[0]
    resid = off - A @ coef
    x2 = float(np.sum(resid ** 2 * w))
    pval = float(chi2d.sf(x2, len(rows) - 2))
    info.update(ttv_p=pval, ttv_max_min=float(np.max(np.abs(resid)) * 1440))
    return ev, info


def run_relaying(cmd, label):
    """Run the fit, passing its progress lines (stderr, "MCMC: ...") to the debug log as
    they come, so a long fit shows how close it is to converging."""
    import threading
    from types import SimpleNamespace
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    err = []

    def pump():
        for line in proc.stderr:
            err.append(line)
            if line.startswith("MCMC:"):
                logging.debug(f"  {label}: {line.strip()}")
    t = threading.Thread(target=pump, daemon=True)
    t.start()
    out = proc.stdout.read()
    proc.wait()
    t.join()
    return SimpleNamespace(returncode=proc.returncode, stdout=out, stderr="".join(err))


def fit_one(r, outdir, steps):
    import find_planets as fp
    import lightkurve as lk
    name, k = r["target"], int(r["candidate"])
    P, T0 = float(r["period_d"]), float(r["t0_btjd"])
    dur = r.get("dur_fit_h")
    D = float(dur if pd.notna(dur) and str(dur) != "" else r["duration_h"]) / 24
    base = dict(target=name, candidate=k)
    star = retry(lambda: stellar_params(name), f"{name}: TIC query")
    if star is None or not (np.isfinite(star["rad"]) and np.isfinite(star["mass"])):
        return dict(base, verdict="no_stellar_params")
    sr = retry(lambda: lk.search_lightcurve(name, author=fp.AUTHORS), f"{name}: MAST search")
    rows, sources = fp.choose_products(sr)
    lcs = [c for c in (fp.clean(lc) for lc in fp.download(sr[rows])) if c is not None]
    ev = local_events(lcs, P, T0, D)
    ttv = {}
    if len(ev) >= 2:
        _, ttv = sector_timing(ev, P, T0, D, float(r["depth_ppm"]) * 1e-6)
    if len(ev) < 2:
        return dict(base, verdict="fit_failed", detail="fewer than 2 transits with data on both sides")
    safe = re.sub(r"[^\w.-]+", "_", name)
    npz = outdir / f"{safe}_c{k}_events.npz"
    np.savez(npz, t=np.concatenate([x[0] for x in ev]), y=np.concatenate([x[1] for x in ev]),
             e=np.concatenate([x[2] for x in ev]), cad=np.concatenate([x[3] for x in ev]))
    rho_tic = star["mass"] / star["rad"] ** 3 * RHO_SUN
    u1, u2 = ld_prior(star["teff"])
    par = dict(P=P, T0=T0, D=D, depth=float(r["depth_ppm"]) * 1e-6, u1=u1, u2=u2, max_steps=steps, rho_star=rho_tic)
    pj = outdir / f"{safe}_c{k}_params.json"
    pj.write_text(json.dumps(par))
    logging.info(f"{name} #{k}: fitting {len(ev)} transits ({sources})")
    res = run_relaying([str(VET_PY), str(HERE / "vetting" / "transit_mcmc.py"), str(npz), str(pj)],
                       f"{name} #{k}")
    lines = [l for l in res.stdout.splitlines() if l.startswith("{")]
    if res.returncode != 0 or not lines:
        return dict(base, verdict="fit_failed", detail=(res.stderr.strip().splitlines() or ["?"])[-1][:200])
    f = json.loads(lines[-1])
    lr, lr_lo, lr_hi = f["log10_rho"]
    # TIC uncertainty in log10 density, from mass and radius errors (default 0.1 dex)
    em = star["e_mass"] / star["mass"] if np.isfinite(star["e_mass"]) and star["e_mass"] > 0 else 0.1
    er = star["e_rad"] / star["rad"] if np.isfinite(star["e_rad"]) and star["e_rad"] > 0 else 0.07
    s_tic = np.hypot(em, 3 * er) / np.log(10)
    diff = np.log10(rho_tic) - lr
    s_fit = (lr_hi - lr) if diff > 0 else (lr - lr_lo)
    z = diff / np.hypot(s_fit, s_tic)
    verdict = "consistent" if abs(z) <= DENSITY_SIGMA else "density_mismatch"
    note = ""
    if verdict == "density_mismatch" and np.isfinite(ttv.get("ttv_p", np.nan)) and ttv["ttv_p"] < TTV_P:
        note = (f"sector timing scatters (offsets up to {ttv['ttv_max_min']:.0f} min, p = {ttv['ttv_p']:.0e}): "
                "if real, timing variations smear the dip and lower the fitted density; for weak "
                "signals the offsets may be noise")
    if verdict == "density_mismatch" and "v_shaped" in str(r.get("flags", "")):
        verdict = "grazing_density_uncertain"
        note = "; ".join(x for x in (f"V-shaped (grazing?): density {z:+.1f} sigma not trusted", note) if x)
    max_t14 = fp.max_duration(f["period"][0], star, f["rp_rs"][0]) * 24
    if f["t14_h"][0] > DUR_MARGIN * max_t14:
        note = "; ".join(x for x in (f"fitted T14 {f['t14_h'][0]:.1f} h vs at most {max_t14:.1f} h for the star "
                                     f"(density {'consistent' if verdict == 'consistent' else 'mismatch'}, {z:+.1f} sigma)", note) if x)
        verdict = "duration_too_long"
    if f.get("converged") is False:
        note = "; ".join(x for x in (note, f"MCMC not converged after {f['steps']} steps "
                                           f"(autocorrelation time {f['autocorr_steps']:.0f})") if x)
    out = dict(base, verdict=verdict, density_z=round(float(z), 2), rho_tic=round(rho_tic, 3),
               n_events=len(ev), steps=f["steps"], autocorr_steps=round(f["autocorr_steps"], 0),
               converged=f.get("converged"),
               sources=sources, events_file=npz.name, ttv_p=ttv.get("ttv_p"),
               ttv_max_min=ttv.get("ttv_max_min"), max_t14_h=round(max_t14, 2), detail=note)
    for key in ("period", "t0", "rp_rs", "b", "rho", "a_rs", "t14_h", "depth_ppm"):
        mid, lo, hi = f[key]
        out[key], out[f"{key}_lo"], out[f"{key}_hi"] = mid, lo, hi
    out["rp_earth"] = out["rp_rs"] * star["rad"] * 109.076
    logging.info(f"  -> {verdict}: rho {out['rho']:.2f} (+{out['rho_hi'] - out['rho']:.2f}/"
                 f"-{out['rho'] - out['rho_lo']:.2f}) vs TIC {rho_tic:.2f} g/cm3 ({z:+.1f} sigma); "
                 f"Rp {out['rp_earth']:.2f} R_earth, b {out['b']:.2f}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", help="enriched.csv or candidates_all.csv")
    ap.add_argument("--target")
    ap.add_argument("--candidate", type=int)
    ap.add_argument("--out")
    ap.add_argument("--steps", type=int, default=30000)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--all", action="store_true", help="every candidate except same_period, not just clean")
    args = ap.parse_args()
    logsetup.setup()
    logging.getLogger("lightkurve").setLevel(logging.WARNING)
    src = Path(args.candidates)
    outdir = Path(args.out) if args.out else src.parent / "transit_fit"
    outdir.mkdir(parents=True, exist_ok=True)
    if not VET_PY.exists():
        sys.exit(f"vetting environment not found ({VET_PY}); see requirements-vet.txt")
    df = pd.read_csv(src)
    if args.target:
        df = df[df.target == args.target]
        if args.candidate:
            df = df[df.candidate == args.candidate]
    elif args.all:
        df = df[df.status != "same_period"]
    else:
        df = df[df.status == "clean"]
    logging.info(f"{len(df)} candidate(s) to fit")
    if df.empty:
        if not args.target:
            (outdir / "transit_fit.csv").unlink(missing_ok=True)   # don't leave an earlier run's results
        return 0

    def safe_fit(r):
        try:
            return fit_one(r, outdir, args.steps)
        except Exception as e:
            logging.exception(f"  {r['target']} #{r['candidate']} failed")
            return dict(target=r["target"], candidate=int(r["candidate"]), verdict="fit_failed",
                        detail=f"{type(e).__name__}: {e}"[:200])

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        rows = list(ex.map(safe_fit, [r for _, r in df.iterrows()]))
    pd.DataFrame(rows).to_csv(outdir / "transit_fit.csv", index=False)
    logging.info(f"wrote {outdir / 'transit_fit.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
