#!/usr/bin/env python3
"""
find_planets.py - unattended BLS transit search over TESS light curves.

Usage:
    python find_planets.py targets/controls.txt [--out results/<name>]
    python find_planets.py --queue results/<name>/queue.db --worker K --out results/<name>/wK
        (worker mode, started by run_parallel.py: claims stars from the run's queue)

For each target in the list (one per line; '#' starts a comment) the script:
  1. downloads all TESS SPOC 2-minute light curves and stitches them,
  2. strips masked/NaN values (avoids the astropy masked-array errors),
  3. iteratively: flattens with known transits masked, trims data near segment
     edges, removes isolated one-cadence dips, removes known transits,
     runs a blind BLS search, vets the best signal, logs it, masks it,
     until the best peak falls below SDE_STOP or SNR_STOP, or MAX_CANDIDATES is reached.
     Signals driven by a few events (glitches) are masked as single events, not logged.

Changes from the earlier version (find_planets_v2.py, tag pre-cleanup):
  - Data are trimmed from both ends of every data segment: at least EDGE_TRIM
    days, more (up to EDGE_TRIM_MAX) while the edge still deviates from the
    segment's level - scattered-light ramps can last a day or more. Ramps and
    poor detrending near gaps made false dips; candidates whose transits mostly
    sit next to a gap are flagged transits_at_data_edges.
  - Single-cadence low outliers (far below the median, neighbours not low) are removed.
  - The glitch test drops the strongest LOO_FRAC of events (at least one), not just
    one, so two glitches that happen to line up at some period are caught.
  - few_transits is flagged below FEW_TRANSITS usable transits, and the odd/even
    test uses a lower threshold for those candidates.
  - The star's radius and mass are looked up in the TIC: candidates bigger than
    MAX_RP_EARTH are flagged too_large_for_planet, and durations longer than the
    star allows are flagged duration_too_long_for_star.
  - The depth is measured in each sector separately (sector_depths_ppm column):
    depth_varies_by_sector and dip_missing_in_sector catch blended eclipsing
    binaries whose depth changes with each pipeline's aperture.
  - Per-sector timing (sector_timing_min column): after removing a linear drift,
    dips that wander by more than TIMING_TOL durations are flagged
    timing_varies_by_sector - typical of evolving starspots, not eclipses.
  - Odd/even mismatches are also flagged double_period_variable, with the doubled
    period in period_x2_d (an eclipsing binary, or two spot groups on a star
    rotating in 2P - interesting either way).
  - Period uniqueness (sde_margin column, "margin@rival period"): the SDE lead over the
    best peak at an unrelated period (not a ratio m/n, m, n <= 5); below
    SDE_MARGIN it is flagged ambiguous_period.
  - Model-shift uniqueness (model_shift column, "ms2;ms3"): folded at the period,
    no other dip or brightening may be about as significant as the candidate,
    allowing for the number of independent phases; otherwise not_unique.
  - The depth at phase 0.5 is measured for every candidate (secondary_ppm); a deep,
    significant one is flagged significant_secondary (a stellar companion).
  - Momentum dumps: per sector, the transits near a dump are counted (dump_near);
    when most are, near_momentum_dumps is flagged (a likely pointing artefact).
  - Stars that log no signal have their downloaded light curves deleted from
    lightkurve's cache afterwards (--keep-lightcurves to keep them); stars with
    any signal keep theirs for the follow-up checks.
  - Logged signals are masked over their measured dip width (mask_width_h), up to
    half the orbit, not just the BLS box, so wide eclipses and spot minima don't
    leave wings that are found again as extra same_period rows.

Outputs (in --out):
  candidates.csv   one row per candidate, with vetting flags
  plots/<target>/  one diagnostic PNG per candidate
  done.txt         finished targets + status (used to resume; errors are retried)
Exit status: 0 when every target has a final status, 3 if any ended in an error
(rerun to retry them), 1 or 128+signal if the script itself failed.
  run.log          progress and error log
"""
import argparse
import csv
import json
import logging
import os
import re
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                      # no display needed on a headless server
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.timeseries import BoxLeastSquares
from astroquery.mast import Catalogs

warnings.filterwarnings("ignore", message=".*oktopus.*")
import lightkurve as lk  # noqa: E402
from mastutil import set_mast_timeout  # noqa: E402
set_mast_timeout()                    # give up on a stuck MAST query after 90 s, not 600

logging.getLogger("lightkurve").setLevel(logging.ERROR)

# ------------------------------------------------------------------ settings
AUTHORS        = ("SPOC", "TESS-SPOC", "QLP")   # per sector, first available wins:
                                # SPOC 2-min target data, else full-frame-image light
                                # curves (200 s from sector 56, 600 s from 27, 1800 s before)
MIN_EXPTIME    = 60             # seconds; skips SPOC 20-s data (huge files, no benefit here)
QLP_STRICT     = True           # drop every QLP cadence with a non-zero quality flag: up to
                                # half a sector can be flagged, and the flagged points carry
                                # most of the outliers and fake deep "transits"
WINDOW_DAYS    = 1.0            # flatten window - must be several x the longest transit
CLIP_SIGMA     = 4.0            # upper-only outlier clip after flattening
BIN_MINUTES    = 10             # bin before searching (5x fewer points, much faster BLS)
P_MIN          = 0.3            # days (ultra-short-period planets go down to ~0.3 d)
P_MAX_CAP      = 50.0           # days; also limited to (data span / 3)
DURATIONS      = np.array([0.5, 0.75, 1, 1.5, 2, 3, 4, 6]) / 24   # days
DUR_COEF       = 0.12           # max trial duration = DUR_COEF * P**(1/3) days
                                # (~2.9 h at 1 d, ~6 h at 10 d): stops absurdly long
                                # boxes at short periods fitting stellar variability
N_CHUNKS       = 12             # period grid is searched in chunks, each with its
                                # own allowed durations
DOWNLOAD_TRIES = 3              # MAST sometimes returns transient 500 errors
OVERSAMPLE     = 3              # period-grid density factor
MAX_GRID       = 400_000        # safety cap on number of trial periods
MAX_CANDIDATES = 6              # per star
MAX_REJECTS    = 8              # glitches + same-period signals handled per star before giving up
GLITCH_HALF    = 0.5            # days masked either side of a rejected single event (minimum)
GLITCH_MERGE   = 1.5            # days: a new glitch this close to an old one widens that mask
ALIAS_TOL      = 0.001          # fractional tolerance for period-ratio alias checks
CONSIST_Z      = 3.0            # per-transit consistency: min z-score of the dip count
CONSIST_MIN_N  = 10             # only apply the consistency test with this many transits
SMALL_N_FRAC   = 0.5            # below CONSIST_MIN_N transits: reject unless MORE than
                                # this fraction dip (noise alone dips in about half)
MAX_JUNK_RUN   = 3              # stop a star after this many non-clean results in a row
                                # (rejections, flagged or same-period signals). Strong
                                # real planets are found first, so this costs little.
SNR_STOP       = 7.0            # stop when the red-noise SNR is below this (cf. Kepler's 7.1)
SDE_STOP       = 9.0            # stop when the best peak's SDE is below this. White noise
                                # alone reaches ~7; calibrate on known planets
MASK_PAD       = 1.5            # mask this multiple of the transit duration
EDGE_TRIM      = 0.5            # days always dropped from each end of every data segment...
EDGE_TRIM_MAX  = 3.0            # ...extended up to this while the edge still looks disturbed
EDGE_BIN       = 0.25           # days; bin width used to judge the edges
EDGE_SIGMA     = 3.0            # an edge bin this many sigma from the segment's level is disturbed
EDGE_GAP       = 0.5            # days; a gap longer than this starts a new segment
LOW_SIGMA      = 5.0            # one-cadence dips this far below the median are removed...
LOW_NEIGHBOUR  = 2.0            # ...unless a neighbouring cadence is this many sigma low too
LOO_FRAC       = 0.25           # glitch test drops this fraction of the strongest events
                                # (at least one); the SNR left must reach SNR_STOP
                                # scaled by sqrt(transits left / transits)
FEW_TRANSITS   = 5              # flag candidates with fewer usable transits than this
ODD_EVEN_SIG   = 3.0            # odd/even depth difference flagged above this many sigma
ODD_EVEN_SIG_FEW = 2.0          # ...or this many, for candidates with few transits
ONE_SIDED_FRAC = 0.5            # flag when more than this fraction of transits lack
                                # data on one side
MIN_BIN_POINTS = 20             # plots: bin so each bin averages about this many points...
MIN_BINS, MAX_BINS = 6, 60      # ...within these limits (bins across +-3 durations)
MAX_RP_EARTH   = 22.4          # 2 Jupiter radii: anything larger is almost always a star
DUR_MARGIN     = 1.5            # flag durations above this multiple of the longest
                                # possible for the catalogued star (allows for catalogue
                                # errors and the coarse DURATIONS grid)
SECTOR_Z       = 5.0            # a sector's depth this many sigma from the mean...
SECTOR_REL     = 0.3            # ...and differing by this fraction is flagged. Pipelines
                                # correct for neighbouring stars differently, so up to
                                # ~15-20% differences are normal
TIMING_Z       = 5.0            # sector dips at least this significant get a timing...
TIMING_TOL     = 0.5            # ...and are flagged if, after removing a linear drift (a
                                # slightly wrong period), they sit more than this many
                                # durations from the prediction. Eclipses and transits keep
                                # time; starspot minima wander as spots evolve
SEC_SIGMA      = 5.0            # a dip at phase 0.5 this significant...
SEC_REL        = 0.1            # ...and deeper than this fraction of the primary...
SEC_ABS_PPM    = 1000           # ...or than this many ppm is a stellar secondary eclipse
DUR_MARGIN_FIT = 1.3            # with a fitted duration (not the coarse BLS grid), flag
                                # durations above this multiple of the longest possible
ASYM_SIGMA     = 3.0            # asymmetric_dip: the two sides differ by this many sigma...
ASYM_FRAC      = 0.25           # ...and by this fraction of the depth
FIT_STARTS     = ((1.0, 0.3),) + tuple((m, f) for m in (0.5, 1.0, 1.5, 2.0, 2.5) for f in (0.3, 0.8))
                                # shape-fit starting points (duration multiple, ingress fraction);
                                # the best fit wins. A single start at the BLS duration can stop
                                # in a short local minimum (TIC 111784356: 3.75 h vs 4.95 h)
SDE_MARGIN     = 1.0            # ambiguous_period: the best peak whose period isn't a simple
                                # ratio (m/n, m, n <= 5) of the candidate's is within this much
                                # SDE of it (TIC 14566446: 10.79 vs 10.51 at an unrelated period)
RIVAL_TOL      = 0.01           # fractional period tolerance for "same or simple ratio"
KEEP_LIGHTCURVES = False        # --keep-lightcurves: don't delete no-signal stars' downloads
MS_STEP        = 0.25           # model-shift test: box depth every this many durations of phase
                                # (at least 100 steps per orbit, so short periods are covered)
OE_TIME_SIGMA  = 3.0            # odd_even_timing: odd and even mid-times differ by this
OE_TIME_FRAC   = 0.15           # many sigma and by this fraction of the duration
DUMP_PAD_H     = 3.0            # a transit within half its duration plus this many hours of a
                                # momentum dump (pointing disturbance) counts as near one
DUMP_FRAC      = 0.5            # flag when more than this fraction of transits (>= 3) are
                                # (hot Jupiters' secondaries are a few hundred ppm at most)

CSV_FIELDS = ["target", "candidate", "period_d", "t0_btjd", "duration_h", "depth_ppm",
              "snr", "snr_white", "snr_without_best_events", "dips", "consistency_z", "sde", "n_transits", "odd_even_sigma", "one_sided_frac", "rp_earth", "max_dur_h", "sector_depths_ppm", "sector_timing_min", "period_x2_d", "mask_width_h", "secondary_ppm", "dump_near", "dur_fit_h", "asym_sigma", "oe_dt_min", "sde_margin", "model_shift", "status", "flags",
              "n_sectors", "sources", "span_d", "plot", "timestamp"]


# ------------------------------------------------------------------ helpers
def arr(q):
    """Any lightkurve/astropy column -> plain float ndarray (masked values -> NaN)."""
    v = getattr(q, "value", q)
    if hasattr(v, "filled"):
        v = v.filled(np.nan)
    return np.asarray(v, dtype=float)


def clean(lc):
    """Rebuild one sector's light curve from plain arrays, normalised to 1,
    dropping masked/NaN points. Some full-frame-image products lack usable
    errors; those get an error estimated from the point-to-point scatter.
    QLP marks bad cadences with its own high quality bits, which lightkurve's
    default mask ignores, so for QLP only unflagged cadences are kept."""
    t, f, e = arr(lc.time), arr(lc.flux), arr(lc.flux_err)
    good = np.isfinite(t) & np.isfinite(f) & (f > 0)
    if QLP_STRICT and lc.meta.get("AUTHOR") == "QLP" and "quality" in lc.colnames:
        good &= arr(lc.quality) == 0
    t, f, e = t[good], f[good], e[good]
    if len(t) < 50:
        return None
    med = np.median(f)
    f, e = f / med, e / med
    if np.mean(np.isfinite(e) & (e > 0)) < 0.9:
        e = np.full(len(f), 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2))
    ok = np.isfinite(e) & (e > 0)
    meta = {k: lc.meta.get(k) for k in ("SECTOR", "AUTHOR")}   # for per-sector checks
    return lk.LightCurve(time=Time(t[ok], format=lc.time.format, scale=lc.time.scale),
                         flux=f[ok], flux_err=e[ok], meta=meta)


def star_params(name, raise_errors=False):
    """Catalogued radius and mass (solar units) from the TIC, or None.
    Needed for the duration and planet-size checks, which are skipped without it.
    raise_errors=True (prefetch.py) raises a failed lookup instead of returning
    None, so "no values in the TIC" and "lookup failed" can be told apart."""
    try:
        m = re.fullmatch(r"TIC\s*(\d+)", name.strip(), flags=re.IGNORECASE)
        if m:
            tab = Catalogs.query_criteria(catalog="TIC", ID=int(m.group(1)))
        else:
            tab = Catalogs.query_object(name, catalog="TIC", radius=0.002 * u.deg)
            tab.sort("dstArcSec")
        if len(tab) == 0:
            return None
        row = tab[0]
        vals = {k: row[k] for k in ("rad", "mass")}
        if any(np.ma.is_masked(v) or not np.isfinite(float(v)) or float(v) <= 0
               for v in vals.values()):
            return None
        return {k: float(v) for k, v in vals.items()}
    except Exception as e:
        if raise_errors:
            raise
        logging.warning(f"  TIC lookup failed ({type(e).__name__}); "
                        f"skipping duration and size checks")
        return None


_manifest = None


def prefetched(name):
    """The star's prefetch.py record, if the run has a manifest (PREFETCH_MANIFEST,
    set by run_parallel.py) with this star fetched and its files still present;
    else None, and the star is searched the usual way."""
    global _manifest
    if _manifest is None:
        _manifest = {}
        path = os.environ.get("PREFETCH_MANIFEST")
        if path and Path(path).exists():
            for line in Path(path).read_text().splitlines():
                try:                                    # a line still being written: skip it
                    r = json.loads(line)
                    _manifest[r["star"]] = r
                except (ValueError, KeyError):
                    continue
    rec = _manifest.get(name)
    if not rec:
        return None
    if rec.get("status") in ("no_data", "no_usable_data"):
        return rec
    if rec.get("status") == "ok" and all(Path(f).is_file() for f in rec.get("files", [])):
        return rec
    return None


def max_duration(P, star, k):
    """Longest possible transit (central, circular orbit) in days, for a body of
    radius ratio k around the catalogued star."""
    a = 215.032 * (star["mass"] * (P / 365.25) ** 2) ** (1 / 3)   # solar radii
    return P / np.pi * np.arcsin(min(1.0, (1 + k) * star["rad"] / a))


def transit_mask(t, planets, pad=MASK_PAD):
    m = np.zeros(len(t), dtype=bool)
    for P, T0, D in planets:
        m |= np.abs((t - T0 + P / 2) % P - P / 2) < 0.5 * D * pad
    return m


def bin_lc(t, y, dy, minutes=BIN_MINUTES):
    w = minutes / 1440
    idx = np.floor((t - t.min()) / w).astype(np.int64)
    _, inv, n = np.unique(idx, return_inverse=True, return_counts=True)
    return (np.bincount(inv, t) / n,
            np.bincount(inv, y) / n,
            np.sqrt(np.bincount(inv, dy ** 2)) / n)


def edge_trim(t, y, gap=EDGE_GAP):
    """True for points kept after trimming segment edges (segments split at gaps
    longer than `gap`). At least EDGE_TRIM days go from each end; walking inward
    in EDGE_BIN bins, trimming continues while bins deviate from the segment's
    level by more than EDGE_SIGMA (its own binned scatter), up to EDGE_TRIM_MAX.
    It stops at two good bins in a row, so a transit inside a good stretch is
    kept. `t` must be sorted; `y` is the flattened flux."""
    keep = np.ones(len(t), dtype=bool)
    for seg in np.split(np.arange(len(t)), np.flatnonzero(np.diff(t) > gap) + 1):
        if len(seg) == 0:
            continue
        ts, ys = t[seg], y[seg]
        span = ts[-1] - ts[0]
        if span < 2 * EDGE_TRIM:
            keep[seg] = False
            continue
        nb = max(1, int(np.ceil(span / EDGE_BIN)))
        idx = np.minimum(((ts - ts[0]) / EDGE_BIN).astype(int), nb - 1)
        cnt = np.bincount(idx, minlength=nb)
        mean = np.bincount(idx, ys, minlength=nb) / np.maximum(cnt, 1)
        centre = (np.arange(nb) + 0.5) * EDGE_BIN
        interior = (cnt > 0) & (centre > EDGE_TRIM_MAX) & (span - centre > EDGE_TRIM_MAX)
        if interior.sum() < 5:                     # short segment: judge against all bins
            interior = cnt > 0
        ref = mean[interior]
        level = np.median(ref)
        sig = 1.4826 * np.median(np.abs(ref - level)) or np.std(ref) or 1e-6
        bad = (cnt == 0) | (np.abs(mean - level) > EDGE_SIGMA * sig)

        def depth_from(order):
            last_bad, run = -1, 0
            for j, b in enumerate(order):
                if (j + 1) * EDGE_BIN > EDGE_TRIM_MAX:
                    return EDGE_TRIM_MAX
                if bad[b]:
                    last_bad, run = j, 0
                else:
                    run += 1
                    if run >= 2 and (j + 1) * EDGE_BIN >= EDGE_TRIM:
                        break
            return max(EDGE_TRIM, (last_bad + 1) * EDGE_BIN)

        lo = depth_from(range(nb))
        hi = depth_from(range(nb - 1, -1, -1))
        keep[seg] = (ts - ts[0] >= lo) & (ts[-1] - ts >= hi)
    return keep


def flatten_sector(lc, protect):
    t_all = arr(lc.time)
    cadence = np.nanmedian(np.diff(t_all))
    window = max(5, int(round(WINDOW_DAYS / cadence)) | 1)  # odd number of cadences
    # flatten() splits the data at gaps longer than break_tolerance cadences. Masked
    # transits look like gaps, so without this a short-period planet's mask chops the
    # light curve into pieces shorter than the window (bad fits, or a crash).
    widest = max((D * MASK_PAD for _, _, D in protect), default=0.0)
    tolerance = max(5, int(np.ceil(widest / cadence)) + 5)
    flat = lc.flatten(window_length=window, mask=transit_mask(t_all, protect),
                      break_tolerance=tolerance)
    t, y, dy = arr(flat.time), arr(flat.flux), arr(flat.flux_err)
    keep = edge_trim(t, y)
    return t[keep], y[keep], dy[keep]


def detrend(lcs, found, protect=None):
    """Flatten each sector (with its own cadence) with transits protected,
    clip, remove `found` transits, bin. `protect` (default: found) lists
    transits kept out of the trend fit."""
    protect = found if protect is None else protect
    parts = [flatten_sector(lc, protect) for lc in lcs]
    t, y, dy = (np.concatenate(p) for p in zip(*parts))
    order = np.argsort(t)
    t, y, dy = t[order], y[order], dy[order]
    ok = np.isfinite(t) & np.isfinite(y) & np.isfinite(dy)
    t, y, dy = t[ok], y[ok], dy[ok]
    med = np.median(y)
    sig = 1.4826 * np.median(np.abs(y - med))
    ok = y < med + CLIP_SIGMA * sig                           # upper outliers
    # One-cadence dips: a transit spans several cadences, so a deep point whose
    # neighbours aren't low is a glitch (a neighbour across a gap counts as not low).
    low = y < med - LOW_NEIGHBOUR * sig
    lone = (y < med - LOW_SIGMA * sig) & ~np.r_[False, low[:-1]] & ~np.r_[low[1:], False]
    ok &= ~lone
    t, y, dy = t[ok], y[ok], dy[ok]
    keep = ~transit_mask(t, found)
    return bin_lc(t[keep], y[keep], dy[keep])


def period_grid(t):
    """Log-spaced periods: keeps transit drift across the data span below
    ~1/OVERSAMPLE of the shortest trial duration at every period."""
    span = t.max() - t.min()
    pmax = min(P_MAX_CAP, span / 3)
    if pmax <= P_MIN * 1.5:
        return None, False
    n = int(OVERSAMPLE * span / DURATIONS.min() * np.log(pmax / P_MIN))
    capped = n > MAX_GRID
    n = int(np.clip(n, 2000, MAX_GRID))
    return np.exp(np.linspace(np.log(P_MIN), np.log(pmax), n)), capped


class Periodogram:
    """Plain container for the combined chunked BLS result."""
    def __init__(self, **kw):
        self.__dict__.update(kw)


def run_bls(bls, periods):
    """BLS in period chunks, each only trying durations physically plausible
    for its shortest period (max duration ~ P**(1/3))."""
    edges = np.exp(np.linspace(np.log(periods[0]), np.log(periods[-1]), N_CHUNKS + 1))
    edges[-1] = np.inf
    keys = ("period", "power", "transit_time", "duration", "depth", "depth_snr")
    parts = {k: [] for k in keys}
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = periods[(periods >= lo) & (periods < hi)]
        if len(sel) == 0:
            continue
        durs = DURATIONS[DURATIONS <= DUR_COEF * sel[0] ** (1 / 3)]
        if len(durs) == 0:
            durs = DURATIONS[:1]
        r = bls.power(sel, durs)
        for k in keys:
            parts[k].append(np.asarray(getattr(r, k), dtype=float))
    return Periodogram(**{k: np.concatenate(v) for k, v in parts.items()})


def local_sde(periods, power, nbins=200):
    """Signal detection efficiency after subtracting the periodogram's running
    median. BLS power climbs towards long periods, so without this, long-period
    noise scores as if it were a signal. (Dividing by local scatter as well was
    tried and rejected: at long periods the power is nearly flat between
    neighbouring trial periods, the local scatter collapses and scores explode.)"""
    nbins = max(10, min(nbins, len(power) // 50))
    lp = np.log(periods)
    edges = np.linspace(lp[0], lp[-1], nbins + 1)
    idx = np.clip(np.digitize(lp, edges) - 1, 0, nbins - 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    med = np.array([np.median(power[idx == b]) for b in range(nbins)])
    base = np.interp(lp, centres, med)
    r = power - base
    return (r - r.mean()) / r.std(), base


def binned(x, y, edges, min_n=1):
    """Bin centres, means and standard errors, for bins with at least min_n points."""
    n, _ = np.histogram(x, edges)
    s, _ = np.histogram(x, edges, weights=y)
    s2, _ = np.histogram(x, edges, weights=y ** 2)
    c = 0.5 * (edges[1:] + edges[:-1])
    ok = n >= max(1, min_n)
    mean = s[ok] / n[ok]
    var = np.maximum(s2[ok] / n[ok] - mean ** 2, 0) * n[ok] / np.maximum(n[ok] - 1, 1)
    return c[ok], mean, np.sqrt(var / n[ok])


def adaptive_edges(lo, hi, n_points):
    """Equal-width bins, as many as the data allow with ~MIN_BIN_POINTS per bin."""
    n = int(np.clip(n_points // MIN_BIN_POINTS, MIN_BINS, MAX_BINS))
    return np.linspace(lo, hi, n + 1)


def dip_width(t, y, P, T0, D):
    """Full width (days) over which the folded flux stays significantly below
    the median, walking outward from mid-transit; at least D. The BLS box can be
    much narrower than a wide eclipse or a spot minimum (durations are capped
    at short periods), and masking only the box leaves wings to be found again."""
    ph = (t - T0 + 0.5 * P) % P - 0.5 * P
    reach = min(0.5 * P, 6 * D)
    w = np.abs(ph) < reach
    if w.sum() < 2 * MIN_BIN_POINTS:
        return D
    n = int(np.clip(w.sum() // MIN_BIN_POINTS, MIN_BINS, max(MIN_BINS, 8 * reach / D)))
    edges = np.linspace(-reach, reach, n + 1)          # bins no narrower than D/4
    bx, by, be = binned(ph[w], y[w], edges, MIN_BIN_POINTS // 2)
    if len(bx) == 0:
        return D
    base = np.median(y)
    j = int(np.argmin(np.abs(bx)))
    # significantly low AND at least 10% of the central depth, so a planet's
    # faint phase curve in very precise data doesn't count as part of the dip
    low = by < base - np.maximum(2 * be, 0.1 * (base - by[j]))
    if not low[j]:
        return D
    a, b = j, j
    while a > 0 and low[a - 1]:
        a -= 1
    while b < len(bx) - 1 and low[b + 1]:
        b += 1
    half = max(abs(bx[a]), abs(bx[b])) + 0.5 * (edges[1] - edges[0])
    return max(D, 2 * half)


def mask_duration(c):
    """Duration to mask for a logged signal: its measured dip width (allowing for
    MASK_PAD), but never masking more than half the orbit."""
    P, T0, D = c["mask"] if c["kind"] == "same_period" else (c["P"], c["T0"], c["D"])
    W = dip_width(c["t"], c["y"], P, T0, D)
    return P, T0, max(D, min(W / MASK_PAD, 0.5 * P / MASK_PAD))


# ------------------------------------------------------------------ vetting
def noise_on_timescale(t, y, P, T0, D):
    """Scatter of the out-of-transit flux averaged over the transit duration.
    Unlike the per-point errors, this includes correlated ("red") noise -
    stellar granulation, spacecraft systematics - which dominates for bright
    stars and makes white-noise SNRs wildly optimistic for long, shallow dips."""
    oot = ~transit_mask(t, [(P, T0, D)])
    tt, yy = t[oot], y[oot]
    idx = np.floor((tt - tt.min()) / D).astype(np.int64)
    _, inv, n = np.unique(idx, return_inverse=True, return_counts=True)
    means = np.bincount(inv, yy) / n
    expected = D / np.median(np.diff(tt))
    means = means[n >= 0.5 * expected]                     # well-filled bins only
    if len(means) < 20:
        return np.nan
    return 1.4826 * np.median(np.abs(means - np.median(means)))


def transit_consistency(t, y, P, T0, D):
    """Count how many individual transits actually show a dip.
    A real planet's dip appears in (nearly) every transit; noise and
    systematics produce dips in only about half. This counting test makes no
    assumption about the noise, so a handful of large glitches can't fake it.
    Returns (number of usable transits, number with a dip, z-score vs 50%)."""
    epoch = np.round((t - T0) / P)
    ph = t - (T0 + epoch * P)
    inn = np.abs(ph) < 0.5 * D
    near = (np.abs(ph) > 0.75 * D) & (np.abs(ph) < 2.5 * D)
    sel = inn | near
    if not sel.any():
        return 0, 0, 0.0
    ep, inv = np.unique(epoch[sel], return_inverse=True)
    ii, yy = inn[sel], y[sel]
    n_in = np.bincount(inv, ii)
    n_out = np.bincount(inv, ~ii)
    s_in = np.bincount(inv, yy * ii)
    s_out = np.bincount(inv, yy * ~ii)
    cadence = np.median(np.diff(t))
    need = max(2, 0.5 * D / cadence)
    ok = (n_in >= need) & (n_out >= need)
    if not ok.any():
        return 0, 0, 0.0
    d = s_out[ok] / n_out[ok] - s_in[ok] / n_in[ok]
    n, pos = int(ok.sum()), int(np.sum(d > 0))
    return n, pos, float((pos - n / 2) / np.sqrt(n / 4))


def trapezoid_drop(dt, depth, t14, tau):
    """Fractional flux drop of a trapezoidal transit; dt from mid-transit."""
    a = np.abs(dt)
    half = t14 / 2
    out = np.where(a <= half - tau, depth, 0.0)
    ramp = (a > half - tau) & (a < half)
    return np.where(ramp, depth * (half - a) / max(tau, 1e-6), out)


def fit_shape(t, y, dy, P, T0, D, depth0):
    """Trapezoid fit to the dip: dict(dtc, P, depth, t14, tau, chi2r), or None.
    The BLS duration comes from a coarse grid, and folding years of data at the
    (capped) grid period smears the dip; so the period and epoch are refined too,
    within the drift the grid allows, and the model is averaged over each bin.
    The fit starts from each of FIT_STARTS and keeps the best."""
    from scipy.optimize import least_squares
    n = np.round((t - T0) / P)
    ph0 = t - T0 - n * P
    w = np.abs(ph0) < min(max(3 * D, 0.25), 0.45 * P)
    if w.sum() < 30 or depth0 <= 0:
        return None
    tt, nn, yy, ee = t[w], n[w], y[w], dy[w]
    span = max(t.max() - t.min(), P)
    dp_max = max(0.5 * D * P / span, 1e-7)
    sub = (np.arange(5) - 2) / 5 * BIN_MINUTES / 1440          # bin-averaged model

    def model(p):
        dtc, dp, depth, t14, f = p
        x = tt - (T0 + dtc) - nn * (P + dp)
        return 1 - np.mean([trapezoid_drop(x + s_, depth, t14, f * t14 / 2) for s_ in sub], axis=0)

    lo = [-D / 2, -dp_max, 0, 0.2 * D, 0.02]
    hi = [D / 2, dp_max, 5 * depth0, min(3 * D, 0.4 * P), 1.0]
    best = None
    for m, f0 in FIT_STARTS:                    # several starts; the lowest chi2 wins
        x0 = [0, 0, depth0, float(np.clip(m * D, lo[3] * 1.01, hi[3] * 0.99)), f0]
        try:
            f = least_squares(lambda p: (yy - model(p)) / ee, x0, bounds=(lo, hi))
        except Exception:
            continue
        if best is None or np.sum(f.fun ** 2) < np.sum(best.fun ** 2):
            best = f
    if best is None:
        return None
    f = best
    dtc, dp, depth, t14, frac = f.x
    return dict(dtc=dtc, P=P + dp, depth=depth, t14=t14, tau=frac * t14 / 2,
                chi2r=np.sum(f.fun ** 2) / max(1, len(yy) - 5))


def shape_tests(t, y, dy, P, T0, fit):
    """Symmetry (sigma of the larger side difference, both conditions applied by
    the caller) and odd/even mid-time difference (minutes, error)."""
    tc, P = T0 + fit["dtc"], fit.get("P", P)
    ph = (t - tc + 0.5 * P) % P - 0.5 * P
    h, d = fit["t14"] / 2, fit["depth"]
    def diff(a_mask, b_mask):
        if a_mask.sum() < 3 or b_mask.sum() < 3:
            return 0.0, np.inf
        sa = np.std(y[a_mask]) / np.sqrt(a_mask.sum())
        sb = np.std(y[b_mask]) / np.sqrt(b_mask.sum())
        return np.mean(y[a_mask]) - np.mean(y[b_mask]), np.hypot(sa, sb)
    # just outside: pre-ingress vs post-egress; inside: first half vs second half
    d_out, e_out = diff((ph > -h - 0.5 * h * 2) & (ph < -h), (ph > h) & (ph < h + 0.5 * h * 2))
    d_in, e_in = diff((ph > -h) & (ph < -0.1 * h), (ph > 0.1 * h) & (ph < h))
    asym = max(abs(d_out) / e_out if (abs(d_out) > ASYM_FRAC * d) else 0,
               abs(d_in) / e_in if (abs(d_in) > ASYM_FRAC * d) else 0)
    # odd/even timing: shift scan with the fitted shape
    epoch = np.round((t - tc) / P).astype(int)
    shifts = np.linspace(-h, h, 81)
    times = []
    for parity in (0, 1):
        sel = (epoch % 2 == parity) & (np.abs(ph) < 3 * h)
        if sel.sum() < 20:
            return asym, np.nan, np.nan
        c2 = np.array([np.sum(((y[sel] - (1 - trapezoid_drop(ph[sel] - s, d, fit["t14"], fit["tau"])))
                               / dy[sel]) ** 2) for s in shifts])
        i = int(np.argmin(c2))
        red = c2[i] / max(1, sel.sum() - 1)
        inside = shifts[c2 <= c2[i] + max(red, 1.0)]
        times.append((shifts[i], (inside.max() - inside.min()) / 2 if len(inside) > 1 else h))
    dt = times[1][0] - times[0][0]
    return asym, dt * 1440, np.hypot(times[0][1], times[1][1]) * 1440


RATIOS = np.array(sorted({m / n for m in range(1, 6) for n in range(1, 6)}))


def related_period(p, P):
    """True where period(s) p are P or a simple ratio m/n of it (m, n <= 5)."""
    r = np.atleast_1d(np.asarray(p, dtype=float)) / P
    return np.min(np.abs(r[:, None] / RATIOS[None, :] - 1), axis=1) < RIVAL_TOL


def rival_peak(periods, z, P, others=()):
    """The strongest periodogram peak (SDE, period) whose period is not P or a
    simple ratio m/n of it (m, n <= 5), nor related to the periods in others
    (signals already found in the star): a competing explanation of the dips.
    A real signal should stand clear of it; noise dips can be strung together
    at several unrelated periods almost equally well."""
    periods = np.asarray(periods)
    related = related_period(periods, P)
    for Po in others:
        related |= related_period(periods, Po)
    if related.all():
        return np.nan, np.nan
    i = int(np.argmax(np.where(related, -np.inf, z)))
    return float(z[i]), float(periods[i])


def model_shift(t, y, P, T0, D):
    """Model-shift uniqueness, after Kepler's Robovetter (Coughlin et al. 2016).
    Fold at P and measure a box depth of duration D at every phase; significances
    come from the scatter of those depths, so correlated noise is included. The
    primary (at T0) is compared with the strongest other dip away from it and from
    the next strongest ("ter"; the strongest, "sec", may be a secondary eclipse)
    and with the strongest brightening ("pos"). FA2 = sqrt(2) erfcinv(D/P) allows
    for the number of independent phases. Returns (ms2, ms3):
      ms2 = sig_pri - sig_ter - FA2   (< 0: another dip is about as significant)
      ms3 = sig_pri - sig_pos - FA2   (< 0: a brightening is about as significant)
    or (nan, nan) when too few phases are covered."""
    from scipy.special import erfcinv
    phd = ((t - T0) / P + 0.5) % 1 * P - 0.5 * P
    o = np.argsort(phd)
    x = np.concatenate([phd[o] - P, phd[o], phd[o] + P])
    cs = np.r_[0, np.cumsum(np.tile(y[o], 3))]
    c = np.arange(-0.5 * P, 0.5 * P, min(MS_STEP * D, P / 100))   # >= 100 phases even for short P
    lo, hi = np.searchsorted(x, c - D / 2), np.searchsorted(x, c + D / 2)
    n = hi - lo
    ok = n >= 0.5 * np.median(n[n > 0]) if (n > 0).any() else n > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        zraw = (np.median(y) - (cs[hi] - cs[lo]) / n) * np.sqrt(n)
    def dist(a, b):                              # phase distance, wrapping around
        d = np.abs(a - b) % P
        return np.minimum(d, P - d)
    away = ok & (dist(c, 0) > 2 * D)
    near = ok & (dist(c, 0) <= 0.25 * D)
    if away.sum() < 20 or not near.any():
        return np.nan, np.nan
    med = np.median(zraw[away])
    unit = 1.4826 * np.median(np.abs(zraw[away] - med))
    if not unit > 0:
        return np.nan, np.nan
    z = (zraw - med) / unit
    pri = np.max(z[near])
    i_sec = np.argmax(np.where(away, z, -np.inf))
    rest = away & (dist(c, c[i_sec]) > 2 * D)
    ter = np.max(z[rest]) if rest.any() else np.nan     # no room for another dip: no ms2
    pos = -np.min(z[away])
    fa2 = np.sqrt(2) * erfcinv(min(D / P, 1.0))
    return float(pri - ter - fa2), float(pri - pos - fa2)


def secondary_depth(t, y, P, T0, D):
    """Depth at phase 0.5 (where a circular binary's second eclipse falls) and
    its error: (depth, err), NaNs without data there."""
    dt = (t - T0) % P - 0.5 * P                    # time from phase 0.5
    inn = np.abs(dt) < 0.5 * D
    near = (np.abs(dt) > 0.75 * D) & (np.abs(dt) < 2.5 * D)
    if inn.sum() < 3 or near.sum() < 6:
        return np.nan, np.nan
    rest = y[~inn]
    sig = 1.4826 * np.median(np.abs(rest - np.median(rest)))
    return (float(np.mean(y[near]) - np.mean(y[inn])),
            float(sig * np.sqrt(1 / inn.sum() + 1 / near.sum())))


def one_sided_fraction(t, P, T0, D):
    """Fraction of transits with in-transit data but too little data just
    before or just after - i.e. sitting against a gap, where detrending is weakest."""
    epoch = np.round((t - T0) / P)
    ph = t - (T0 + epoch * P)
    inn = np.abs(ph) < 0.5 * D
    before = (ph < -0.75 * D) & (ph > -2.5 * D)
    after = (ph > 0.75 * D) & (ph < 2.5 * D)
    sel = inn | before | after
    if not inn.any():
        return 0.0
    _, inv = np.unique(epoch[sel], return_inverse=True)
    n_in = np.bincount(inv, inn[sel])
    n_b = np.bincount(inv, before[sel])
    n_a = np.bincount(inv, after[sel])
    need = max(2, 0.5 * D / np.median(np.diff(t)))
    has = n_in >= need
    if not has.any():
        return 0.0
    return float(np.mean((n_b[has] < need) | (n_a[has] < need)))


def sector_depths(t, y, P, T0, D, sectors):
    """Depth measured separately in each sector: (label, depth, error) for sectors
    with enough data in and beside the transits. A planet's depth is the same in
    every sector; a blended eclipsing binary's changes with each pipeline's
    aperture, and mixed data can hide a dip that is missing from some sectors."""
    ph = (t - T0 + 0.5 * P) % P - 0.5 * P
    inn = np.abs(ph) < 0.5 * D
    near = (np.abs(ph) > 0.75 * D) & (np.abs(ph) < 2.5 * D)
    need = max(3, 0.5 * D / np.median(np.diff(t)))
    out = []
    for label, lo, hi in sectors:
        s = (t >= lo) & (t <= hi)
        i, o = s & inn, s & near
        if i.sum() < need or o.sum() < 2 * need:
            continue
        rest = y[s & ~inn]
        sig = 1.4826 * np.median(np.abs(rest - np.median(rest)))
        d = np.mean(y[o]) - np.mean(y[i])
        out.append((label, float(d), float(sig * np.sqrt(1 / i.sum() + 1 / o.sum()))))
    return out


def sector_timing(t, y, P, T0, D, W, sectors):
    """Dip-centre offset per sector: (label, offset in days, mid-time), for
    sectors whose dip is at least TIMING_Z sigma. The offset is the centroid of
    the flux deficit within +-W/2 (W = measured dip width), re-centred a few times
    so a shifted dip isn't truncated by the window."""
    ph = (t - T0 + 0.5 * P) % P - 0.5 * P
    half = 0.5 * max(W, D)
    need = max(3, 0.5 * D / np.median(np.diff(t)))
    out = []
    for label, lo, hi in sectors:
        s = (t >= lo) & (t <= hi)
        p, ys = ph[s], y[s]
        inn, outside = np.abs(p) < 0.5 * D, np.abs(p) >= half
        if inn.sum() < need or outside.sum() < 2 * need:
            continue
        base = np.median(ys[outside])
        sig = 1.4826 * np.median(np.abs(ys[outside] - base))
        if base - np.mean(ys[inn]) < TIMING_Z * sig / np.sqrt(inn.sum()):
            continue
        c = 0.0
        for _ in range(3):
            q = (p - c + 0.5 * P) % P - 0.5 * P
            win = np.abs(q) < half
            wts = base - ys[win]
            if wts.sum() <= 0:
                break
            c += np.sum(q[win] * wts) / wts.sum()
        out.append((label, float(c), 0.5 * (lo + hi)))
    return out


def timing_varies(series, D):
    """True if, after removing a linear drift (a slightly wrong period), any
    sector's dip sits more than TIMING_TOL durations from the prediction."""
    if len(series) < 3:                 # 2 points always fit a line exactly
        return False
    tm = np.array([x[2] for x in series])
    off = np.array([x[1] for x in series])
    resid = off - np.polyval(np.polyfit(tm, off, 1), tm)
    return bool(np.max(np.abs(resid)) > TIMING_TOL * D)


def red_snr(depth, n_tr, sigma_d):
    return depth / sigma_d * np.sqrt(n_tr) if np.isfinite(sigma_d) and sigma_d > 0 else np.nan


def vet(P, D, depth, stats, found, n_c, one_sided, star, per_sector):
    """n_c: transits usable by the dip-count test; one_sided: fraction of
    transits against a gap; star: catalogued radius/mass or None;
    per_sector: output of sector_depths()."""
    flags = []
    if star and depth > 0:
        k = np.sqrt(depth)
        if k * star["rad"] * 109.1 > MAX_RP_EARTH:
            flags.append("too_large_for_planet")
        if D > DUR_MARGIN * max_duration(P, star, k):
            flags.append("duration_too_long_for_star")
    if len(per_sector) >= 2:
        d = np.array([x[1] for x in per_sector])
        e = np.array([x[2] for x in per_sector])
        dw = np.sum(d / e ** 2) / np.sum(1 / e ** 2)          # weighted mean depth
        if np.any((np.abs(d - dw) > SECTOR_Z * e) & (np.abs(d - dw) > SECTOR_REL * abs(dw))):
            flags.append("depth_varies_by_sector")
        if np.any((dw > 3 * e) & (d < e)):                   # should be seen, isn't
            flags.append("dip_missing_in_sector")
    for Pf, _, _ in found:                       # alias (m/n ratio) of a found planet?
        r = P / Pf
        if any(abs(r - m / n) < ALIAS_TOL * (m / n)
               for m in range(1, 6) for n in range(1, 6) if np.gcd(m, n) == 1):
            flags.append(f"alias_of_{Pf:.4f}")
    n_tr = int(np.sum(np.asarray(stats["per_transit_count"]) > 0))
    few = n_c < FEW_TRANSITS            # e.g. 4/4 dips says little
    if few:
        flags.append("few_transits")
    if one_sided > ONE_SIDED_FRAC:
        flags.append("transits_at_data_edges")
    (do, eo), (de, ee) = stats["depth_odd"], stats["depth_even"]
    oe_sig = abs(do - de) / np.hypot(eo, ee) if np.hypot(eo, ee) > 0 else 0.0
    rel = abs(do - de) / max(abs(do), abs(de), 1e-12)
    # With few transits the odd/even errors are large, so use a lower bar.
    limit = ODD_EVEN_SIG_FEW if few else ODD_EVEN_SIG
    if oe_sig > limit and rel > 0.1:    # at very high SNR, tiny differences are "significant"
        flags.append("odd_even_mismatch")
        # The true period is probably twice this: an eclipsing binary's unequal
        # eclipses, or two spot groups on a star rotating in 2P. Worth a look either way.
        flags.append("double_period_variable")
    ll = np.asarray(stats["per_transit_log_likelihood"])
    pos = ll[ll > 0]
    if n_tr >= 2 and pos.sum() > 0 and pos.max() / pos.sum() > 0.5:
        flags.append("single_event_dominated")
    if stats["harmonic_delta_log_likelihood"] > 0:
        flags.append("sinusoidal_variability")
    if D / P > 0.1 and not star:        # fallback; with stellar parameters the
        flags.append("long_duration")       # duration_too_long_for_star test applies
                                            # (fixed 10% wrongly flags real P < 1 d planets)
    if depth <= 0:
        flags.append("non_positive_depth")
    return flags, n_tr, oe_sig


# ------------------------------------------------------------------ plotting
def save_plot(path, name, k, c):
    res, t, y = c["res"], c["t"], c["y"]
    P, T0, D = c["P"], c["T0"], c["D"]
    snr, sde, flags = c["snr"], c["sde"], c["flags"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 11))

    ax = axes[0]
    ax.plot(res.period, c["z"], lw=0.5, color="k")
    ax.axvline(P, color="r", alpha=0.3, lw=3)
    ax.axhline(SDE_STOP, color="tab:blue", ls="--", lw=1, label=f"SDE_STOP = {SDE_STOP}")
    ax.set_xscale("log")
    ax.set_xlabel("Period [d]")
    ax.set_ylabel("SDE")
    ax.set_title("Periodogram, relative to local baseline")
    ax.legend(loc="upper left")

    ph = (t - T0 + 0.5 * P) % P - 0.5 * P
    w = np.abs(ph) < 3 * D
    min_n = MIN_BIN_POINTS // 2           # drop bins left nearly empty by gaps

    ax = axes[1]
    ax.scatter(ph[w] * 24, y[w], s=2, color="0.5", alpha=0.5)
    edges = adaptive_edges(-3 * D, 3 * D, w.sum())
    bx, by, be = binned(ph[w], y[w], edges, min_n)
    ax.errorbar(bx * 24, by, be, color="r", lw=2, capsize=2)
    ax.set_xlabel("Hours from mid-transit")
    ax.set_ylabel("Normalised flux")
    ax.set_title(f"Folded ({len(edges) - 1} bins)")

    ax = axes[2]
    epoch = np.round((t - T0) / P).astype(int)
    parities = ((0, "tab:blue", "even"), (1, "tab:orange", "odd"))
    # shared edges, sized for the sparser parity, so the two compare like for like
    fewest = min(np.sum(w & (epoch % 2 == p)) for p, _, _ in parities)
    edges = adaptive_edges(-3 * D, 3 * D, fewest)
    for parity, colour, label in parities:
        sel = w & (epoch % 2 == parity)
        if sel.any():
            bx, by, be = binned(ph[sel], y[sel], edges, min_n)
            ax.errorbar(bx * 24, by, be, color=colour, lw=2, capsize=2, label=label)
    ax.legend()
    ax.set_xlabel("Hours from mid-transit")
    ax.set_ylabel("Normalised flux")
    ax.set_title(f"Odd vs even transits (should match for a planet; {len(edges) - 1} bins)")

    fig.suptitle(f"{name}  #{k}:  P = {P:.5f} d,  dur = {D * 24:.2f} h,  "
                 f"SNR = {snr:.1f},  SDE = {sde:.1f}\n"
                 f"flags: {', '.join(flags) if flags else 'none'}")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ------------------------------------------------------------------ search
def dump_times(lc):
    """Momentum-dump times (quality bit 32, reaction-wheel desaturation) from a
    light curve's file, re-read without the quality mask that removed them."""
    try:
        raw = lk.read(lc.meta["FILENAME"], quality_bitmask="none")
        q = arr(raw.quality).astype(np.int64)
        t = arr(raw.time)[(q & 32) != 0]
        return t[np.r_[True, np.diff(t) > 0.1]] if len(t) else t
    except Exception:
        return np.array([])


def dump_proximity(t, P, T0, D, dumps, sectors):
    """Per sector: transits with data, and how many lie within D/2 + DUMP_PAD_H of
    a momentum dump. Returns (near, total, text)."""
    dumps = np.asarray(dumps, float)
    ph = (t - T0 + 0.5 * P) % P - 0.5 * P
    epochs = np.unique(np.round((t[np.abs(ph) < 0.5 * D] - T0) / P))
    tc = T0 + epochs * P
    near_total, total, parts = 0, 0, []
    for label, lo, hi in sectors:
        sel = tc[(tc >= lo) & (tc <= hi)]
        if len(sel) == 0:
            continue
        near = 0
        if len(dumps):
            near = int(np.sum(np.min(np.abs(dumps[None, :] - sel[:, None]), axis=1)
                              < 0.5 * D + DUMP_PAD_H / 24))
        near_total, total = near_total + near, total + len(sel)
        parts.append(f"{label}={near}/{len(sel)}")
    return near_total, total, ";".join(parts)


def search_once(lc, found, excluded, k, star=None, sectors=(), dumps=()):
    """One search pass. Returns None (stop), {"kind": "reject", ...} for a
    single-event glitch, or {"kind": "candidate", ...}."""
    masked = found + excluded
    t, y, dy = detrend(lc, masked)
    if len(t) < 100:
        logging.info("  too few points left after masking; stopping")
        return None
    periods, capped = period_grid(t)
    if periods is None:
        logging.info("  data span too short to search; stopping")
        return None
    if capped and not found and not excluded:          # first pass only: once per star
        # Worst-case transit-time misalignment across the span (half a grid step),
        # against the design target; the same at every period, as the grid is log-spaced.
        span = t.max() - t.min()
        smear = 0.5 * np.log(periods[1] / periods[0]) * span * 1440
        target = 0.5 * DURATIONS.min() / OVERSAMPLE * 1440
        logging.warning(f"  period grid capped at {MAX_GRID} points (span {span:.0f} d): "
                        f"transits misalign by up to {smear:.0f} min vs {target:.0f} min "
                        f"target; the shortest durations are slightly less sensitive")

    res = run_bls(BoxLeastSquares(t, y, dy), periods)
    z, _ = local_sde(res.period, res.power)
    i = int(np.argmax(z))
    sde = float(z[i])
    P, T0, D = float(res.period[i]), float(res.transit_time[i]), float(res.duration[i])
    logging.info(f"  #{k}: P = {P:.5f} d, dur = {D * 24:.2f} h, "
                 f"depth = {res.depth[i] * 1e6:.0f} ppm, SNR = {res.depth_snr[i]:.1f}, "
                 f"SDE = {sde:.1f}")
    if sde < SDE_STOP:
        logging.info(f"  SDE below {SDE_STOP}; stopping")
        return None

    # Re-flatten with this candidate's transits also kept out of the trend fit,
    # then re-measure: the first pass erodes the depth (the "shoulders" problem).
    t, y, dy = detrend(lc, masked, protect=masked + [(P, T0, D)])
    # A dip at the edge of a data segment can vanish here (edge trimming); with no
    # in-transit data left compute_stats fails. It was an edge artefact: mask it.
    n_in = int(np.sum(np.abs((t - T0 + 0.5 * P) % P - 0.5 * P) < 0.5 * D))
    bls = BoxLeastSquares(t, y, dy)
    try:
        stats = bls.compute_stats(P, D, T0) if n_in >= 3 else None
    except ValueError:
        stats = None
    if stats is None:
        logging.info(f"     only {n_in} in-transit point(s) left after re-detrending (dip at a "
                     f"data edge); masking it")
        return dict(kind="inconsistent", mask=(P, T0, D))
    depth, depth_err = (float(v) for v in stats["depth"])
    snr_white = depth / depth_err if depth_err > 0 else 0.0
    n_tr = int(np.sum(np.asarray(stats["per_transit_count"]) > 0))
    sigma_d = noise_on_timescale(t, y, P, T0, D)
    snr = red_snr(depth, n_tr, sigma_d)
    # The red-noise estimate can only make the SNR more cautious, never less.
    snr = min(snr, snr_white) if np.isfinite(snr) else snr_white
    logging.info(f"     re-measured: depth = {depth * 1e6:.0f} ppm, SNR = {snr:.1f} "
                 f"(white-noise SNR {snr_white:.1f})")
    if snr < SNR_STOP:
        logging.info(f"  SNR below {SNR_STOP}; stopping")
        return None

    # Leave-some-out: a planet's signal survives losing its strongest transits
    # (scaled for the transits left); a few glitches that line up do not.
    ll = np.asarray(stats["per_transit_log_likelihood"], dtype=float)
    times = np.asarray(stats["transit_times"], dtype=float)
    n_drop = min(n_tr, max(1, int(round(LOO_FRAC * n_tr))))
    top = np.argsort(ll)[::-1][:n_drop]
    top = top[ll[top] > 0] if np.any(ll[top] > 0) else top[:1]   # only events that dip
    tcs = times[top]
    n_left = n_tr - len(tcs)
    keep = np.all(np.abs(t[:, None] - tcs[None, :]) > D * MASK_PAD / 2, axis=1)
    try:
        d2 = float(BoxLeastSquares(t[keep], y[keep], dy[keep])
                   .compute_stats(P, D, T0)["depth"][0])
        snr_loo = red_snr(d2, n_left, sigma_d) if np.isfinite(d2) and n_left > 0 else 0.0
        alt = snr * d2 / depth * np.sqrt(n_left / n_tr) if depth > 0 else 0.0
        snr_loo = min(snr_loo, alt) if np.isfinite(snr_loo) else alt
    except ValueError:              # no in-transit data left: it was only those events
        snr_loo = 0.0
    need = SNR_STOP * np.sqrt(max(n_left, 0) / n_tr) if n_tr > 0 else SNR_STOP
    if not snr_loo >= need:
        when = ", ".join(f"{x:.2f}" for x in tcs)
        logging.info(f"     driven by {len(tcs)} event(s) at {when} (SNR without them = "
                     f"{snr_loo:.1f}, needed {need:.1f}); masking them and searching again")
        return dict(kind="reject", tcs=tcs, half=max(GLITCH_HALF, D))

    n_c, dips, cz = transit_consistency(t, y, P, T0, D)
    common = dict(P=P, T0=T0, D=D, depth=depth, snr=snr, snr_white=snr_white,
                  snr_loo=snr_loo, sde=sde, n_tr=n_tr, res=res, z=z, t=t, y=y,
                  dips=f"{dips}/{n_c}", cz=cz)

    # Related to a planet already found: a secondary eclipse (found at P, or at
    # P/2, P/3 because the primaries are masked) or a phase-curve harmonic
    # (any small-integer ratio). Logged, masked, and not counted as a planet.
    for Pf, T0f, _ in found:
        what = None
        for n in (1, 2, 3):
            if abs(n * P / Pf - 1) < ALIAS_TOL:
                phases = [((T0 + j * P - T0f) / Pf) % 1.0 for j in range(n)]
                j = int(np.argmin([abs(p - 0.5) for p in phases]))
                if abs(phases[j] - 0.5) < 0.05:
                    what = "secondary_eclipse"
                    common.update(P=Pf, T0=T0 + j * P)
                else:
                    what = "phase_curve_residual"
                break
        if what is None and any(abs(P / Pf - m / n) < ALIAS_TOL * (m / n)
                                for m in range(1, 6) for n in range(1, 6)
                                if np.gcd(m, n) == 1):
            what = "harmonic"
        if what:
            logging.info(f"     related to the {Pf:.5f} d planet: {what}")
            return dict(kind="same_period", flags=[f"related_to_{Pf:.4f}", what],
                        oe_sig=0.0, one_sided=0.0, per_sector=[], timing={}, star=None,
                        sec=(np.nan, np.nan), dump="", dur_fit=np.nan, asym=np.nan, oe_dt=np.nan,
                        ms=(np.nan, np.nan), rival=(np.nan, np.nan), mask=(P, T0, D), **common)

    if n_c < 3:
        logging.info(f"     only {n_c} transit(s) with data; rejecting and masking it")
        return dict(kind="inconsistent", mask=(P, T0, D))
    if n_c >= CONSIST_MIN_N and cz < CONSIST_Z:
        logging.info(f"     only {dips} of {n_c} individual transits show a dip "
                     f"(z = {cz:.1f}); rejecting as noise and masking it")
        return dict(kind="inconsistent", mask=(P, T0, D))
    if n_c < CONSIST_MIN_N and dips <= SMALL_N_FRAC * n_c:
        logging.info(f"     only {dips} of {n_c} individual transits show a dip; "
                     f"rejecting as noise and masking it")
        return dict(kind="inconsistent", mask=(P, T0, D))

    one_sided = one_sided_fraction(t, P, T0, D)
    per_sector = sector_depths(t, y, P, T0, D, sectors)
    n_near, n_tot, dump_text = dump_proximity(t, P, T0, D, dumps, sectors)
    flags, n_tr, oe_sig = vet(P, D, depth, stats, found, n_c, one_sided, star, per_sector)
    rival_sde, rival_p = rival_peak(res.period, z, P, [f[0] for f in found])
    if np.isfinite(rival_sde) and sde - rival_sde < SDE_MARGIN:
        flags.append("ambiguous_period")             # an unrelated period does about as well
    sec, sec_err = secondary_depth(t, y, P, T0, D)
    ms2, ms3 = model_shift(t, y, P, T0, D)
    if (np.isfinite(ms2) and ms2 < 0) or (np.isfinite(ms3) and ms3 < 0):
        flags.append("not_unique")                  # the fold has other dips (or bumps) this strong
    # measured shape: duration (replacing the coarse BLS value in the stellar test),
    # symmetry, odd/even timing
    fit = fit_shape(t, y, dy, P, T0, D, depth)
    dur_fit, asym, oe_dt, oe_err = np.nan, np.nan, np.nan, np.nan
    if fit:
        dur_fit = fit["t14"]
        if star and fit["depth"] > 0:
            mx = max_duration(P, star, np.sqrt(fit["depth"]))
            flags[:] = [f for f in flags if f != "duration_too_long_for_star"]
            if dur_fit > DUR_MARGIN_FIT * mx:
                flags.append("duration_too_long_for_star")
        asym, oe_dt, oe_err = shape_tests(t, y, dy, P, T0, fit)
        if asym > ASYM_SIGMA:
            flags.append("asymmetric_dip")
        if np.isfinite(oe_dt) and abs(oe_dt) > OE_TIME_SIGMA * oe_err and \
                abs(oe_dt) > OE_TIME_FRAC * dur_fit * 1440:
            flags.append("odd_even_timing")
    if n_tot >= 3 and n_near >= 3 and n_near / n_tot > DUMP_FRAC:
        flags.append("near_momentum_dumps")         # may be a pointing artefact
    if np.isfinite(sec) and sec > SEC_SIGMA * sec_err and \
            (sec > SEC_REL * depth or sec > SEC_ABS_PPM * 1e-6):
        flags.append("significant_secondary")       # a stellar companion's second eclipse
    # Timing per sector. For a double-period variable, time each minimum at 2P:
    # folded at P they overlap, and a fixed one hides a wandering one.
    W = dip_width(t, y, P, T0, D)
    if "double_period_variable" in flags:
        timing = {"A": sector_timing(t, y, 2 * P, T0, D, W, sectors),
                  "B": sector_timing(t, y, 2 * P, T0 + P, D, W, sectors)}
    else:
        timing = {"": sector_timing(t, y, P, T0, D, W, sectors)}
    if any(timing_varies(v, D) for v in timing.values()):
        flags.append("timing_varies_by_sector")
    return dict(kind="candidate", P=P, T0=T0, D=D, depth=depth, snr=snr,
                snr_white=snr_white, snr_loo=snr_loo, sde=sde,
                dips=f"{dips}/{n_c}", cz=cz, flags=flags, n_tr=n_tr, oe_sig=oe_sig,
                one_sided=one_sided, per_sector=per_sector, timing=timing, star=star,
                sec=(sec, sec_err), dump=dump_text, dur_fit=dur_fit, asym=asym, oe_dt=oe_dt,
                ms=(ms2, ms3), rival=(rival_sde, rival_p),
                res=res, z=z, t=t, y=y)


def add_glitch(windows, tc, half):
    """Mask a one-off event. Events close to an existing mask widen it, so an
    extended disturbance isn't rediscovered piece by piece."""
    for i, (c0, h0) in enumerate(windows):
        if abs(tc - c0) < GLITCH_MERGE + h0:
            lo, hi = min(c0 - h0, tc - half), max(c0 + h0, tc + half)
            windows[i] = ((lo + hi) / 2, (hi - lo) / 2)
            return
    windows.append((tc, half))


def search_star(lc, name, n_sectors, outdir, writer, fh, sources="", star=None, sectors=(),
                dumps=()):
    """Iterative BLS search on an already-cleaned light curve."""
    safe = re.sub(r"[^\w.-]+", "_", name)
    plotdir = outdir / "plots" / safe
    found = []          # planet candidates: (P, T0, D)
    glitches = []       # one-off events: (centre, half-width)
    same = []           # same-period signals, masked like transits: (P, T0, D)
    n_rejects = 0
    junk_run = 0        # consecutive non-clean results
    k = 0               # rows logged
    pending = []        # (row, c, plot): written when the star is finished, so that
                        # ambiguous_period can be dropped when the rival period turns
                        # out to be another signal found in this star (a second planet)

    while len(found) < MAX_CANDIDATES and n_rejects < MAX_REJECTS:
        if junk_run >= MAX_JUNK_RUN:
            logging.info(f"  {junk_run} non-clean results in a row; stopping")
            break
        # glitch windows are expressed as a single "transit" with a huge period
        excluded = same + [(1e6, c, 2 * h / MASK_PAD) for c, h in glitches]
        try:
            c = search_once(lc, found, excluded, k + 1, star, sectors, dumps)
        except Exception:
            # Keep the candidates already logged; don't fail (and later retry) the star.
            logging.exception(f"  search pass for #{k + 1} failed; stopping this star")
            write_rows(pending, name, writer, fh)
            return f"partial_{len(found)}_candidates_then_error"
        if c is None:
            break
        if c["kind"] == "reject":
            for tc in c["tcs"]:
                add_glitch(glitches, tc, c["half"])
            n_rejects += 1
            junk_run += 1
            continue
        if c["kind"] == "inconsistent":
            same.append(c["mask"])
            n_rejects += 1
            junk_run += 1
            continue

        k += 1
        plotdir.mkdir(parents=True, exist_ok=True)
        plot = plotdir / f"cand{k}_P{c['P']:.4f}.png"
        save_plot(plot, name, k, c)

        if c["kind"] == "same_period":
            status = "same_period"
        else:
            status = "flagged" if c["flags"] else "clean"
        # ambiguous_period alone doesn't count as junk: the rival may be a second planet
        serious = [f for f in c["flags"] if f != "ambiguous_period"]
        junk_run = 0 if status != "same_period" and not serious else junk_run + 1
        mask = mask_duration(c)
        if mask[2] > 1.05 * c["D"]:
            logging.info(f"     dip is wider than the {c['D'] * 24:.2f} h box; masking "
                         f"{mask[2] * MASK_PAD * 24:.2f} h")
        pending.append(({
            "target": name, "candidate": k, "period_d": f"{c['P']:.6f}",
            "t0_btjd": f"{c['T0']:.5f}", "duration_h": f"{c['D'] * 24:.2f}",
            "depth_ppm": f"{c['depth'] * 1e6:.0f}", "snr": f"{c['snr']:.2f}",
            "snr_white": f"{c['snr_white']:.2f}",
            "snr_without_best_events": f"{c['snr_loo']:.2f}",
            "dips": c["dips"], "consistency_z": f"{c['cz']:.1f}",
            "sde": f"{c['sde']:.2f}", "n_transits": c["n_tr"],
            "odd_even_sigma": f"{c['oe_sig']:.2f}",
            "one_sided_frac": f"{c['one_sided']:.2f}",
            "rp_earth": (f"{np.sqrt(max(c['depth'], 0)) * c['star']['rad'] * 109.1:.1f}"
                         if c["star"] else ""),
            "max_dur_h": (f"{max_duration(c['P'], c['star'], np.sqrt(max(c['depth'], 0))) * 24:.2f}"
                          if c["star"] else ""),
            "sector_depths_ppm": ";".join(f"{x[0]}={x[1] * 1e6:.0f}+-{x[2] * 1e6:.0f}"
                                          for x in c["per_sector"]),
            "sector_timing_min": " | ".join(
                (f"{k}: " if k else "") + ";".join(f"{lab}={o * 1440:+.0f}" for lab, o, _ in v)
                for k, v in c["timing"].items() if v),
            "period_x2_d": (f"{2 * c['P']:.6f}"
                            if "double_period_variable" in c["flags"] else ""),
            "mask_width_h": f"{mask[2] * MASK_PAD * 24:.2f}",
            "dump_near": c["dump"],
            "dur_fit_h": f"{c['dur_fit'] * 24:.2f}" if np.isfinite(c["dur_fit"]) else "",
            "asym_sigma": f"{c['asym']:.1f}" if np.isfinite(c["asym"]) else "",
            "oe_dt_min": f"{c['oe_dt']:+.0f}" if np.isfinite(c["oe_dt"]) else "",
            "sde_margin": (f"{c['sde'] - c['rival'][0]:.2f}@{c['rival'][1]:.4f}"
                           if np.isfinite(c["rival"][0]) else ""),
            "model_shift": (";".join(f"{v:+.1f}" if np.isfinite(v) else "" for v in c["ms"])
                            if np.isfinite(c["ms"]).any() else ""),
            "secondary_ppm": (f"{c['sec'][0] * 1e6:.0f}+-{c['sec'][1] * 1e6:.0f}"
                              if np.isfinite(c["sec"][0]) else ""),
            "status": status, "flags": ";".join(c["flags"]),
            "n_sectors": n_sectors, "sources": sources, "span_d": f"{c['t'].max() - c['t'].min():.1f}",
            "plot": str(plot.relative_to(outdir)),
            "timestamp": datetime.now().isoformat(timespec="seconds"),
        }, c, plot))
        if c["kind"] == "same_period":
            same.append(mask)
            n_rejects += 1
        else:
            found.append(mask)

    write_rows(pending, name, writer, fh)
    if glitches:
        logging.info(f"  masked {len(glitches)} glitch window(s)")
    return f"ok_{len(found)}_candidates"


def write_rows(pending, name, writer, fh):
    """Write a star's rows. A candidate whose rival period (sde_margin) matches
    another signal logged for this star is not ambiguous: it has a neighbour,
    e.g. a second planet, so ambiguous_period is dropped (and the plot redrawn)."""
    for row, c, plot in pending:
        if "ambiguous_period" in c["flags"] and np.isfinite(c["rival"][1]):
            others = [o["P"] for _, o, _ in pending if o is not c]
            if any(related_period(c["rival"][1], Po)[0] for Po in others):
                c["flags"].remove("ambiguous_period")
                row["flags"] = ";".join(c["flags"])
                if row["status"] == "flagged" and not c["flags"]:
                    row["status"] = "clean"
                logging.info(f"     #{row['candidate']}: rival period {c['rival'][1]:.4f} d is another "
                             f"signal in this star; ambiguous_period dropped")
                save_plot(plot, name, int(row["candidate"]), c)
        writer.writerow(row)
    fh.flush()


def download(sr):
    attempt, repaired = 1, set()
    while True:
        try:
            return sr.download_all()
        except Exception as e:
            # lightkurve reuses any cached file that exists; one left half-written by
            # an interrupted download (a worker killed mid-star) fails to read every
            # time. Delete it and download it again, at once (not counted as a try).
            bad = corrupt_cache_file(e)
            if bad and bad not in repaired:
                repaired.add(bad)
                logging.warning(f"  cached file is corrupt (interrupted download?): {bad.name}; "
                                f"deleting it and downloading again")
                bad.unlink()
                continue
            if attempt >= DOWNLOAD_TRIES:
                raise
            logging.warning(f"  download attempt {attempt} failed ({type(e).__name__}); "
                            f"retrying in 60 s")
            attempt += 1
            time.sleep(60)


def corrupt_cache_file(exc):
    """The cached file a lightkurve read error names as possibly corrupt, if it is
    inside lightkurve's cache folder; else None."""
    msg = str(exc)
    if "corrupt" not in msg:
        return None
    cache = Path(lk.config.get_cache_dir()).resolve()
    for m in re.finditer(r"(\S+\.fits(?:\.gz)?)", msg):
        f = Path(m.group(1)).resolve()
        if f.is_file() and cache in f.parents:
            return f
    return None


def choose_products(sr):
    """One light curve per sector, by AUTHORS priority, then shortest cadence."""
    tab = sr.table
    sector = np.asarray(tab["sequence_number"], dtype=int)
    author = np.asarray(tab["author"]).astype(str)
    exptime = np.asarray(tab["exptime"], dtype=float)
    best = {}
    for i in range(len(tab)):
        if author[i] not in AUTHORS or exptime[i] < MIN_EXPTIME:
            continue
        rank = (AUTHORS.index(author[i]), exptime[i])
        if sector[i] not in best or rank < best[sector[i]][0]:
            best[sector[i]] = (rank, i)
    picks = sorted(best.items())
    rows = [i for _, (_, i) in picks]
    summary = {}
    for s, (_, i) in picks:
        summary.setdefault(f"{author[i]}/{int(exptime[i])}s", []).append(str(s))
    label = "; ".join(f"{k}: {','.join(v)}" for k, v in summary.items())
    return rows, label


def process_target(name, outdir, writer, fh):
    rec = prefetched(name)
    if rec:                                             # fetched ahead: no MAST contact
        if rec["status"] != "ok":
            return rec["status"]
        lcc = lk.LightCurveCollection([lk.read(f) for f in rec["files"]])
        sources = rec["sources"]
        logging.info(f"  {len(lcc)} sector(s) prefetched {rec['fetched'][:10]} ({sources})")
    else:
        for attempt in range(1, DOWNLOAD_TRIES + 1):    # MAST searches time out too
            try:
                sr = lk.search_lightcurve(name, author=AUTHORS)
                break
            except Exception as e:
                if attempt == DOWNLOAD_TRIES:
                    raise
                logging.warning(f"  search attempt {attempt} failed ({type(e).__name__}); retrying in 60 s")
                time.sleep(60)
        if len(sr) == 0:
            return "no_data"
        rows, sources = choose_products(sr)
        if not rows:
            return "no_usable_data"
        logging.info(f"  {len(rows)} sector(s) chosen ({sources}); downloading")
        lcc = download(sr[rows])
        if lcc is None or len(lcc) == 0:
            return "error: download_failed"              # usually transient: retried
    lcs = [c for c in (clean(lc) for lc in lcc) if c is not None]
    if not lcs:
        if not KEEP_LIGHTCURVES:
            drop_downloads(lcc)
        return "no_good_data"
    logging.info(f"  {sum(len(c) for c in lcs)} good cadences from {len(lcs)} sector(s)")
    star = rec.get("star_params") if rec else star_params(name)
    if star:
        logging.info(f"  TIC: R = {star['rad']:.2f} R_sun, M = {star['mass']:.2f} M_sun")
    sectors = [(f"{lc.meta.get('SECTOR')}{'(QLP)' if lc.meta.get('AUTHOR') == 'QLP' else ''}",
                float(arr(lc.time).min()), float(arr(lc.time).max())) for lc in lcs]
    dumps = np.concatenate([dump_times(lc) for lc in lcc]) if len(lcc) else np.array([])
    status = search_star(lcs, name, len(lcs), outdir, writer, fh, sources, star, sectors, dumps)
    if status == "ok_0_candidates" and not KEEP_LIGHTCURVES:
        drop_downloads(lcc)
    return status


def drop_downloads(lcc):
    """Delete a star's downloaded light curves from lightkurve's cache (and folders
    left empty). Used for stars that logged no signal, about 90% of them: nothing
    later needs their data, and a re-search simply downloads it again (~6 MB a star).
    Stars with any signal keep theirs for the follow-up checks and later review."""
    cache = Path(lk.config.get_cache_dir()).resolve()
    freed = 0
    for lc in lcc:
        f = Path(str(lc.meta.get("FILENAME", ""))).resolve()
        if not f.is_file() or cache not in f.parents:      # only ever inside the cache
            continue
        freed += f.stat().st_size
        f.unlink()
        d = f.parent
        while d != cache and d.is_dir() and not any(d.iterdir()):
            d.rmdir()
            d = d.parent
    if freed:
        logging.info(f"  no signal: removed its light curves from the cache ({freed / 1e6:.1f} MB)")


# ------------------------------------------------------------------ main
def read_targets(path):
    out = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def read_done(path):
    """Targets already finished. Errors are not counted, so they are retried."""
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2 and not parts[1].startswith("error"):
                done.add(parts[0])
    return done


def drop_rows(csv_path, names):
    """Remove rows for these targets from a candidates file. A star that failed
    partway may have written rows before the error; it is about to be searched
    again, so they go (else it appears twice in the results)."""
    if not csv_path.exists() or not names:
        return
    with open(csv_path, newline="") as f:
        rows = list(csv.DictReader(f))
    keep = [r for r in rows if r.get("target") not in set(names)]
    if len(keep) < len(rows):
        tmp = csv_path.with_suffix(".tmp")
        with open(tmp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
            w.writeheader()
            w.writerows(keep)
        os.replace(tmp, csv_path)
        logging.info(f"removed {len(rows) - len(keep)} row(s) left by an earlier attempt")


def search_one(name, outdir, label=""):
    """Search one star, append its rows to candidates.csv and its outcome to
    done.txt; returns the status. An exception becomes an "error: ..." status, so
    one bad star can't stop the worker."""
    csv_path = outdir / "candidates.csv"
    logging.info(f"{label}{name}")
    start = time.time()
    new_csv = not csv_path.exists() or csv_path.stat().st_size == 0
    with open(csv_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        if new_csv:
            writer.writeheader()
        try:
            status = process_target(name, outdir, writer, fh)
        except Exception as e:
            logging.exception(f"  failed: {name}")
            status = f"error: {type(e).__name__}: {e}".replace("\n", " ")
    logging.info(f"  -> {status} ({(time.time() - start) / 60:.1f} min)")
    with open(outdir / "done.txt", "a") as d:
        d.write(f"{name}\t{status}\t{datetime.now().isoformat(timespec='seconds')}\n")
    return status


def run_queue(queue_path, worker, outdir):
    """Worker mode (run_parallel.py): claim the best-ranked star that is due from
    the run's queue, search it, record the outcome, repeat until nothing is due
    or <run folder>/STOP exists."""
    from workqueue import Queue
    q = Queue(queue_path)
    stop = Path(queue_path).parent / "STOP"
    commit = os.environ.get("RUN_COMMIT")
    n = 0
    while not stop.exists():
        name = q.claim(worker, os.getpid(), commit)
        if name is None:
            break
        n += 1
        drop_rows(outdir / "candidates.csv", [name])
        status = search_one(name, outdir, f"[worker {worker} #{n}] ")
        q.finish(name, status)
    logging.info(f"worker {worker}: {n} star(s) searched; "
                 f"{'STOP file found' if stop.exists() else 'nothing left that is due'}; exiting")
    q.close()
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("targets", nargs="?", help="text file, one target per line (not with --queue)")
    ap.add_argument("--out", default="results/single", help="output folder (default: results/single)")
    ap.add_argument("--queue", help="worker mode: claim stars from this queue.db (run_parallel.py)")
    ap.add_argument("--worker", type=int, default=0, help="worker number (with --queue)")
    ap.add_argument("--keep-lightcurves", action="store_true",
                    help="keep downloaded light curves of stars with no signal (default: delete them)")
    args = ap.parse_args()
    if bool(args.targets) == bool(args.queue):
        ap.error("give a target list or --queue, not both")
    global KEEP_LIGHTCURVES
    KEEP_LIGHTCURVES = args.keep_lightcurves

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(outdir / "run.log"),
                                  logging.StreamHandler()])
    if args.queue:
        return run_queue(args.queue, args.worker, outdir)

    targets = read_targets(args.targets)
    done_path = outdir / "done.txt"
    done = read_done(done_path)
    todo = [t for t in targets if t not in done]
    logging.info(f"{len(targets)} targets, {len(targets) - len(todo)} already done, "
                 f"{len(todo)} to process")
    drop_rows(outdir / "candidates.csv", todo)
    for i, name in enumerate(todo, 1):
        search_one(name, outdir, f"[{i}/{len(todo)}] ")

    # exit status: 3 if any of this list's targets ended in an error (they are
    # retried on a rerun), 0 when every target has a final status
    latest = {}
    for line in done_path.read_text().splitlines() if done_path.exists() else []:
        parts = line.split("\t")
        if len(parts) >= 2:
            latest.setdefault(parts[0], []).append(parts[1])
    failed = [t for t in targets if not any(not st.startswith("error") for st in latest.get(t, []))]
    if failed:
        logging.warning(f"run complete; {len(failed)} target(s) not finished: {', '.join(failed)}")
        return 3
    logging.info("run complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
