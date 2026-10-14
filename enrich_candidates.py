#!/usr/bin/env python3
"""
enrich_candidates.py - add catalogue information to find_planets candidates.

Usage:
    python enrich_candidates.py results/<name>/candidates_all.csv [--out FILE] [--all]

For each star in the candidates file (queried once per star, cached as JSON):
  - ExoFOP-TESS: TOIs, community TOIs, follow-up imaging/spectroscopy/time series,
    known stellar companions
  - SIMBAD at the target: main identifier, object type, spectral type
  - Known variables within VAR_RADIUS_ARCMIN: VSX, Gaia DR3 eclipsing binaries
    and classifications, ZTF periodic variables
  - TIC neighbours within TIC_RADIUS_ARCMIN, with TESS magnitudes
  - Gaia DR3 binarity indicators for the target: RUWE, non-single-star solution,
    radial-velocity variability (Katz et al. 2023 criterion), corrected BP/RP
    flux excess C* (Riello et al. 2021), double-peaked scans

For each candidate row:
  - period_matches: known variables whose period matches P, 2P, P/2 or a
    3:2-family alias (3P/2, 2P/3, 3P, P/3), or a one-day alias of any of these
    (a ground survey's catalogued period can be one). A star
    near the target with the same period is almost certainly the real source
    (contact binaries often come out at half their period).
  - Flags: nearby_variable_period_match (another star shares the period - a
    blend), known_planet_period_match (the target is a catalogued planet host),
    target_known_variable_period_match, target_catalogued_variable (any
    catalogued variable at the target, with or without a period), has_toi,
    has_ctoi, has_spoc_tce (SPOC's own pipeline detected this period, or twice or
    half it, as a TCE: a known signal even if never promoted to a TOI), simbad_<type>,
    gaia_ruwe_high, gaia_non_single_star, gaia_rv_variable, gaia_bp_rp_excess,
    gaia_multi_peak.
  - blend_suspects: neighbours that could produce the observed depth at all.
    The amplitude a neighbour would need is
        observed depth x (target light in aperture + neighbour light leaking in)
                       / neighbour light leaking in,
    using a simple TESS image model (PSF_* settings). Neighbours needing less
    than SUSPECT_MAX_AMP (a deep eclipse) are listed. This scales the search
    with brightness: a bright neighbour can contaminate from much further away.
    Shallow signals nearly always have some suspects (a faint neighbour with a
    deep eclipse), so this is not a flag: it lists which stars a pixel-level
    check should look at, easiest first (min_blend_amp_pct).

Writes <input dir>/enriched.csv (all original columns plus the new ones).
By default same_period rows are skipped (they repeat a logged signal); --all keeps them.
"""
import argparse
import json
import tempfile
import xml.etree.ElementTree as ET
import logging
import re
import time
import warnings
from functools import lru_cache
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord

warnings.filterwarnings("ignore")
from astroquery.mast import Catalogs, Observations   # noqa: E402
from astroquery.simbad import Simbad   # noqa: E402
from mastutil import set_mast_timeout  # noqa: E402
set_mast_timeout()                    # give up on a stuck MAST query after 90 s, not 600
from astroquery.vizier import Vizier   # noqa: E402

# ------------------------------------------------------------------ settings
VAR_RADIUS_ARCMIN = 3.0         # known-variable search radius (~8.5 TESS pixels)
TIC_RADIUS_ARCMIN = 5.0         # neighbours considered as blend sources
SIMBAD_RADIUS_ARCSEC = 10.0
PERIOD_TOL = 0.005              # relative tolerance for a period match
HARMONICS = {1.0: "P", 2.0: "2P", 0.5: "P/2",      # contact binaries often show at P/2
             1.5: "3P/2", 2 / 3: "2P/3", 3.0: "3P", 1 / 3: "P/3"}   # and 3:2 aliases
TARGET_ARCSEC = 5.0             # a catalogue entry this close is the target itself
SIDEREAL_DAY_F = 1.00274         # cycles per day: ground surveys' daily alias spacing
PIXEL_ARCSEC = 21.0             # TESS pixel scale
APERTURE_PX = 2.0               # radius of a typical photometric aperture
PSF_SIGMA_PX = 1.5              # core of the star image (generous: this is a screen)
PSF_HALO = (0.1, 4.0)           # fraction of light, sigma (px) in faint broad wings
SUSPECT_MAX_AMP = 0.5           # neighbours needing less than this amplitude are suspects
MAX_SUSPECTS = 3                # listed per candidate
QUERY_TRIES = 3
GAIA_RUWE_MAX = 1.4             # above: astrometric wobble, often an unresolved companion
GAIA_MULTI_PEAK = 10            # % of scans showing a double image
GAIA_EXCESS_SIGMA = 3.0         # corrected BP/RP flux excess C* (Riello et al. 2021)

VARIABLE_CATALOGUES = {         # VizieR: label -> (catalogue, name column, period function)
    "VSX": ("B/vsx/vsx", "Name", lambda r: r["Period"]),
    "Gaia-EB": ("I/358/veb", "Source", lambda r: 1 / r["Freq"] if r["Freq"] > 0 else np.nan),
    "Gaia-class": ("I/358/vclassre", "Source", lambda r: np.nan),
    "ZTF": ("J/ApJS/249/18", "ID", lambda r: r["Per"]),
    # TESS Ten Thousand Catalog (Kostov et al. 2025, 2025ApJS..279...50K): vetted TESS FFI
    # eclipsing binaries, new (table 3) and known (table 4). Most aren't in VSX yet
    "T10k-new": ("J/ApJS/279/50/table3", "TIC", lambda r: r["Per"], "EB"),
    "T10k-known": ("J/ApJS/279/50/table4", "TIC", lambda r: r["Per-TESS"], "EB"),
}
LATER_CATALOGUES = ("T10k-new", "T10k-known")   # added 2026-10-14: queried for older caches too
EXOFOP_URL = "https://exofop.ipac.caltech.edu/tess/target.php"
TCE_PERIOD_TOL = 0.003          # relative tolerance for matching a SPOC TCE's period


# ------------------------------------------------------------------ helpers
def retry(what, fn):
    for attempt in range(1, QUERY_TRIES + 1):
        try:
            return fn()
        except Exception as e:
            if attempt == QUERY_TRIES:
                logging.warning(f"  {what} failed ({type(e).__name__}: {e}); skipped")
                return None
            time.sleep(20 * attempt)


def num(x):
    """Plain float, NaN for masked/missing values."""
    try:
        return float(np.ma.filled(x, np.nan)) if np.ma.is_masked(x) else float(x)
    except (TypeError, ValueError):
        return np.nan


def resolve_tic(name):
    """TIC row for 'TIC <n>' or any name SIMBAD/MAST can resolve."""
    m = re.fullmatch(r"TIC\s*(\d+)", name.strip(), flags=re.IGNORECASE)
    if m:
        tab = Catalogs.query_criteria(catalog="TIC", ID=int(m.group(1)))
    else:
        tab = Catalogs.query_object(name, catalog="TIC", radius=0.002 * u.deg)
        tab.sort("dstArcSec")
    return tab[0] if len(tab) else None


def exofop_summary(tic):
    d = requests.get(EXOFOP_URL, params={"id": str(tic), "json": ""}, timeout=60).json()
    return {
        "tois": [t.get("toi", "") + (f" ({t['disp_tfop'] or t['disp_tess']})"
                                     if t.get("disp_tfop") or t.get("disp_tess") else "")
                 for t in d.get("tois", [])],
        "ctois": [t.get("ctoi", "") for t in d.get("ctois", [])],
        "n_imaging": len(d.get("imaging", [])),
        "n_spectroscopy": len(d.get("spectroscopy", [])),
        "n_time_series": len(d.get("time_series", [])),
        "companions": [f"{c.get('sep', '')}\" dmag {c.get('dmag', '')}"
                       for c in d.get("stellar_companions", [])],
    }


def spoc_tces(tic):
    """Threshold-crossing events from SPOC's data-validation reports on MAST:
    [{sectors, period, epoch_btjd, mes, eb}]. Reads two reports' XML (~0.3 MB
    each): the widest sector range and the latest, which together cover what
    SPOC has searched. A star with no 2-min SPOC data has none."""
    obs = Observations.query_criteria(obs_collection="TESS", target_name=str(tic),
                                      dataproduct_type="timeseries")
    if len(obs) == 0:
        return []
    prods = Observations.get_product_list(obs)
    xml = prods[[str(f).endswith("_dvr.xml") for f in prods["productFilename"]]]
    if len(xml) == 0:
        return []

    def rng(f):
        m = re.search(r"-s(\d{4})-s(\d{4})-", f)
        return (int(m[1]), int(m[2])) if m else (0, 0)
    names = [str(f) for f in xml["productFilename"]]
    widest = max(names, key=lambda f: (rng(f)[1] - rng(f)[0], rng(f)[1], f))
    latest = max(names, key=lambda f: (rng(f)[1], rng(f)[1] - rng(f)[0], f))
    want = sorted({widest, latest})
    out = []
    with tempfile.TemporaryDirectory() as tmp:
        Observations.download_products(xml[[f in want for f in names]], download_dir=tmp,
                                       cache=False)
        for path in Path(tmp).rglob("*_dvr.xml"):
            lo, hi = rng(path.name)
            for c in ET.parse(path).getroot().iter("{http://www.nasa.gov/2018/TESS/DV}planetCandidate"):
                out.append({"sectors": f"s{lo}-s{hi}", "period": float(c.get("orbitalPeriodInDays")),
                            "epoch_btjd": float(c.get("epochTjd")),
                            "mes": float(c.get("maxMultipleEventSigma")),
                            "eb": c.get("suspectedEclipsingBinary") == "true"})
    return out


def tce_matches(P, tces):
    """SPOC TCEs at the candidate's period, twice or half it."""
    return [t for t in tces if any(abs(t["period"] / (P * m) - 1) < TCE_PERIOD_TOL
                                   for m in (1.0, 2.0, 0.5))]


def simbad_summary(coord):
    s = Simbad()
    s.add_votable_fields("otype", "sp_type")
    tab = s.query_region(coord, radius=SIMBAD_RADIUS_ARCSEC * u.arcsec)
    if tab is None or len(tab) == 0:
        return {}
    sep = SkyCoord(tab["ra"], tab["dec"], unit="deg").separation(coord).arcsec
    # prefer the star over planet entries at the same position ("pi. Men b")
    sep = np.where(np.char.startswith(np.asarray(tab["otype"], str), "Pl"), sep + 1e3, sep)
    r = tab[int(np.argmin(sep))]
    sep = np.where(sep > 1e3, sep - 1e3, sep)
    return {"main_id": str(r["main_id"]), "otype": str(r["otype"]),
            "sp_type": str(r["sp_type"]), "sep": round(float(np.min(sep)), 1)}


def variables_near(coord, labels=None):
    v = Vizier(columns=["**", "+_r"], row_limit=-1)
    out = []
    for label, (cat, name_col, period, *default_type) in VARIABLE_CATALOGUES.items():
        if labels is not None and label not in labels:
            continue
        res = retry(f"{label} query", lambda: v.query_region(
            coord, radius=VAR_RADIUS_ARCMIN * u.arcmin, catalog=cat))
        if not res:
            continue
        # a catalogue can hold several tables (ZTF's has more than one): use the one with
        # the columns we need, not simply the first returned
        tab = None
        for t in res.values():
            if name_col in t.colnames and len(t):
                try:
                    period(t[0])
                except KeyError:
                    continue
                tab = t
                break
        if tab is None:
            continue
        ra_col = next(c for c in ("RAJ2000", "RA_ICRS", "RAdeg") if c in tab.colnames)
        dec_col = ra_col.replace("RA", "DE")
        seps = SkyCoord(np.asarray(tab[ra_col], float), np.asarray(tab[dec_col], float),
                        unit="deg").separation(coord).arcsec
        for r, sep in zip(tab, seps):
            out.append({"cat": label, "name": ("TIC " if name_col == "TIC" else "") + str(r[name_col]),
                        "sep": round(float(sep), 1),
                        "type": str(r["Type"]) if "Type" in r.colnames
                        else str(r["Class"]) if "Class" in r.colnames
                        else default_type[0] if default_type else "",
                        "period": num(period(r))})
    return out


def tic_neighbours(coord, tic):
    tab = Catalogs.query_region(coord, radius=TIC_RADIUS_ARCMIN * u.arcmin, catalog="TIC")
    out = []
    for r in tab:
        if str(r["ID"]) == str(tic) or not np.isfinite(num(r["Tmag"])):
            continue
        out.append({"tic": str(r["ID"]), "tmag": round(num(r["Tmag"]), 2),
                    "sep": round(num(r["dstArcSec"]), 1)})
    return out


GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap"   # ESA Gaia archive


def gaia_query(adql):
    """Run an ADQL query on the Gaia archive (synchronous; up to 2000 rows) and
    return an astropy Table. pyvo, not astroquery.gaia, which contacts the
    archive on import to print its news."""
    import pyvo
    return pyvo.dal.TAPService(GAIA_TAP).run_sync(adql).to_table()


def gaia_id_of(row):
    """The TIC's Gaia source id as text (a 19-digit id doesn't survive a float)."""
    v = row["GAIA"]
    v = "" if np.ma.is_masked(v) else str(v).strip()
    return v if v.isdigit() else None


def gaia_binarity(gaia_id, coord):
    """Gaia DR3 binarity indicators for the target (by the TIC's Gaia id, else
    the nearest source within 2")."""
    cols = ("source_id, phot_g_mean_mag, bp_rp, ruwe, ipd_frac_multi_peak, non_single_star, "
            "rv_nb_transits, rv_chisq_pvalue, rv_renormalised_gof, rv_amplitude_robust, "
            "phot_bp_rp_excess_factor")
    where = (f"source_id = {gaia_id}" if gaia_id else
             f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {coord.ra.deg}, {coord.dec.deg}, "
             f"{2 / 3600}))=1")
    tab = gaia_query(f"SELECT {cols} FROM gaiadr3.gaia_source WHERE {where}")
    if len(tab) == 0:
        return {}
    r = tab[0]
    g = {c.lower(): num(r[c]) for c in tab.colnames if c.lower() != "source_id"}
    g["source_id"] = str(r[tab.colnames[0]])
    x, G, C = g.get("bp_rp"), g.get("phot_g_mean_mag"), g.get("phot_bp_rp_excess_factor")
    if all(np.isfinite(v) for v in (x, G, C)):           # corrected excess, Riello+2021 eq. 6
        f = (1.154360 + 0.033772 * x + 0.032277 * x ** 2 if x < 0.5 else
             1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3 if x < 4 else
             1.057572 + 0.140537 * x)
        g["c_star"] = C - f
        g["c_star_sigma"] = 0.0059898 + 8.817481e-12 * G ** 7.618399
    return {k: (None if isinstance(v, float) and not np.isfinite(v) else v) for k, v in g.items()}


def gaia_flags(g):
    """Flags from gaia_binarity()'s values."""
    if not g:
        return []
    v = lambda k: g.get(k) if g.get(k) is not None else np.nan
    flags = []
    if v("ruwe") > GAIA_RUWE_MAX:
        flags.append("gaia_ruwe_high")
    if v("non_single_star") > 0:
        flags.append("gaia_non_single_star")
    if v("rv_nb_transits") >= 10 and v("rv_chisq_pvalue") <= 0.01 and v("rv_renormalised_gof") > 4:
        flags.append("gaia_rv_variable")                  # Katz et al. 2023
    if abs(v("c_star")) > GAIA_EXCESS_SIGMA * v("c_star_sigma"):
        flags.append("gaia_bp_rp_excess")
    if v("ipd_frac_multi_peak") > GAIA_MULTI_PEAK:
        flags.append("gaia_multi_peak")
    return flags


def fetch_star(name, cache_dir):
    """Everything catalogue-side for one star, cached as JSON."""
    path = cache_dir / (re.sub(r"[^\w.-]+", "_", name) + ".json")
    if path.exists():
        info = json.loads(path.read_text())
        if "gaia" not in info:                           # cache from before the Gaia checks
            coord = SkyCoord(info["ra"] * u.deg, info["dec"] * u.deg)
            row = retry(f"TIC lookup for {name}", lambda: resolve_tic(name))
            info["gaia"] = retry("Gaia DR3", lambda: gaia_binarity(
                gaia_id_of(row) if row is not None else None, coord)) or {}
            path.write_text(json.dumps(info, indent=1))
        if not set(LATER_CATALOGUES) <= set(info.get("catalogues_checked", [])):
            coord = SkyCoord(info["ra"] * u.deg, info["dec"] * u.deg)   # cache from before them
            info["variables"] = [x for x in info.get("variables", []) if x["cat"] not in LATER_CATALOGUES] \
                + variables_near(coord, LATER_CATALOGUES)
            info["catalogues_checked"] = list(VARIABLE_CATALOGUES)
            path.write_text(json.dumps(info, indent=1))
        if "spoc_tces" not in info:                      # cache from before the TCE check
            tces = retry("SPOC TCEs", lambda: spoc_tces(info["tic"]))
            if tces is not None:
                info["spoc_tces"] = tces
                path.write_text(json.dumps(info, indent=1))
        return info
    row = retry(f"TIC lookup for {name}", lambda: resolve_tic(name))
    if row is None:
        return None
    tic = str(row["ID"])
    coord = SkyCoord(num(row["ra"]) * u.deg, num(row["dec"]) * u.deg)
    info = {"tic": tic, "ra": coord.ra.deg, "dec": coord.dec.deg, "tmag": num(row["Tmag"]),
            "exofop": retry("ExoFOP", lambda: exofop_summary(tic)) or {},
            "simbad": retry("SIMBAD", lambda: simbad_summary(coord)) or {},
            "variables": variables_near(coord), "catalogues_checked": list(VARIABLE_CATALOGUES),
            "neighbours": retry("TIC neighbours", lambda: tic_neighbours(coord, tic)) or []}
    info["gaia"] = retry("Gaia DR3", lambda: gaia_binarity(gaia_id_of(row), coord)) or {}
    tces = retry("SPOC TCEs", lambda: spoc_tces(tic))
    if tces is not None:                    # a failed query is retried next time
        info["spoc_tces"] = tces
    path.write_text(json.dumps(info, indent=1))
    return info


@lru_cache(maxsize=None)
def light_in_aperture(offset_px):
    """Fraction of a star's light inside a circular aperture of radius
    APERTURE_PX whose centre is offset_px from the star (core + halo model)."""
    r, h = APERTURE_PX, 0.05
    g = np.arange(-r, r + h, h)
    X, Y = np.meshgrid(g, g)
    inside = X ** 2 + Y ** 2 <= r ** 2
    d2 = (X - offset_px) ** 2 + Y ** 2
    total = 0.0
    for weight, s in ((1 - PSF_HALO[0], PSF_SIGMA_PX), PSF_HALO):
        total += weight * np.sum(np.exp(-d2[inside] / (2 * s * s))) * h * h / (2 * np.pi * s * s)
    return total


def period_matches(P, variables):
    """Known variables with a period matching a HARMONICS multiple of P, as (kind, text).
    Ground-based surveys observe about once a night, so a catalogued period can
    be a one-sidereal-day alias of the true one (frequency |f +- 1.0027| per day);
    those match too. kind: 'nearby' (another star - probably the real source),
    'planet' (the target is a catalogued transiting-planet host, VSX type EP) or
    'target'."""
    out = []
    for v in variables:
        if not np.isfinite(v["period"]) or v["period"] <= 0:
            continue
        fv = 1 / v["period"]
        for h, label in HARMONICS.items():
            fc = 1 / (h * P)                              # candidate frequency at this harmonic
            how = None
            if abs(v["period"] / (h * P) - 1) < PERIOD_TOL:
                how = label
            else:
                for fa in (abs(fc - SIDEREAL_DAY_F), fc + SIDEREAL_DAY_F):
                    # same frequency precision as the direct match: the alias shifts a
                    # frequency, it doesn't make it any more precise
                    if fa > 0 and abs(fv - fa) < PERIOD_TOL * fc:
                        how = f"1-day alias of {label}"
            if how is None:
                continue
            if v["sep"] >= TARGET_ARCSEC:
                kind = "nearby"
            elif re.search(r"\bEP\b", v["type"]):
                kind = "planet"
            else:
                kind = "target"
            out.append((kind, f"{v['cat']} {v['name']} ({v['type'] or 'EB'}, "
                              f"P={v['period']:.5f} d = {how}, {v['sep']:.0f}\")"))
            break
    return out


def blend_suspects(depth, tmag, neighbours):
    """Neighbours that could produce `depth`, easiest first: (amplitude needed, text)."""
    if not (depth > 0 and np.isfinite(tmag)):
        return []
    own = light_in_aperture(0.0)
    out = []
    for n in neighbours:
        if n["sep"] < 2:                          # the target itself / duplicates
            continue
        L = 10 ** (-0.4 * (n["tmag"] - tmag))     # neighbour / target brightness
        leak = L * light_in_aperture(round(n["sep"] / PIXEL_ARCSEC, 2))
        if leak <= 0:
            continue
        need = depth * (own + leak) / leak
        if need < SUSPECT_MAX_AMP:
            out.append((need, f"TIC {n['tic']} (T={n['tmag']:.1f}, {n['sep']:.0f}\", "
                              f"needs {need * 100:.0f}%)"))
    return sorted(out)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", help="candidates.csv or candidates_all.csv")
    ap.add_argument("--out", help="output CSV (default: enriched.csv beside the input)")
    ap.add_argument("--all", action="store_true", help="also enrich same_period rows")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    src = Path(args.candidates)
    df = pd.read_csv(src)
    if not args.all:
        df = df[df["status"] != "same_period"].copy()
    cache_dir = src.parent / "enrich_cache"
    cache_dir.mkdir(exist_ok=True)

    stars = {}
    for i, name in enumerate(df["target"].unique(), 1):
        logging.info(f"[{i}/{df['target'].nunique()}] {name}")
        stars[name] = fetch_star(name, cache_dir)

    cols = {k: [] for k in ("tic_id", "simbad_id", "simbad_type", "simbad_sptype",
                            "exofop_tois", "exofop_ctois", "exofop_followup",
                            "exofop_companions", "spoc_tces", "gaia_ruwe", "gaia_rv_chisq_pvalue",
                            "gaia_c_star_sigmas", "target_variable", "n_variables_nearby",
                            "period_matches",
                            "n_blend_suspects", "min_blend_amp_pct", "blend_suspects",
                            "enrich_flags")}
    for _, r in df.iterrows():
        s = stars.get(r["target"]) or {}
        ex, sb = s.get("exofop", {}), s.get("simbad", {})
        matches = period_matches(float(r["period_d"]), s.get("variables", []))
        suspects = blend_suspects(float(r["depth_ppm"]) * 1e-6, s.get("tmag", np.nan),
                                  s.get("neighbours", []))
        flags = []
        kinds = {k for k, _ in matches}
        for kind, flag in (("nearby", "nearby_variable_period_match"),
                           ("planet", "known_planet_period_match"),
                           ("target", "target_known_variable_period_match")):
            if kind in kinds:
                flags.append(flag)
        # a catalogued variable at the target's own position, period or not (Gaia's
        # classifications often have none): the dip may be the star's own variability
        own = [v for v in s.get("variables", [])
               if v["sep"] < TARGET_ARCSEC and not re.search(r"\bEP\b", v["type"])]
        # Gaia's SOLAR_LIKE class (low-level activity, common in planet hosts) is
        # listed but doesn't raise the flag
        if any(v["type"] != "SOLAR_LIKE" for v in own) and \
                "target_known_variable_period_match" not in flags:
            flags.append("target_catalogued_variable")
        cols["target_variable"].append("; ".join(
            f"{v['cat']} {v['type'] or 'variable'}" + (f" P={v['period']:.5f} d" if np.isfinite(v["period"]) else "")
            for v in own))
        if ex.get("tois"):
            flags.append("has_toi")
        if ex.get("ctois"):
            flags.append("has_ctoi")
        tces = tce_matches(float(r["period_d"]), s.get("spoc_tces", []))
        if tces:
            flags.append("has_spoc_tce")
        cols["spoc_tces"].append("; ".join(
            f"P={t['period']:.5f} d MES {t['mes']:.1f} ({t['sectors']}{', EB?' if t['eb'] else ''})"
            for t in tces))
        if "*" in sb.get("otype", "") and sb.get("otype") not in ("*", "PM*", "**"):
            flags.append(f"simbad_{sb['otype']}")
        g = s.get("gaia", {})
        flags += gaia_flags(g)
        cols["gaia_ruwe"].append(f"{g['ruwe']:.2f}" if g.get("ruwe") is not None else "")
        cols["gaia_rv_chisq_pvalue"].append(f"{g['rv_chisq_pvalue']:.1e}"
                                            if g.get("rv_chisq_pvalue") is not None else "")
        cols["gaia_c_star_sigmas"].append(f"{g['c_star'] / g['c_star_sigma']:.1f}"
                                          if g.get("c_star") is not None else "")
        cols["tic_id"].append(s.get("tic", ""))
        cols["simbad_id"].append(sb.get("main_id", ""))
        cols["simbad_type"].append(sb.get("otype", ""))
        cols["simbad_sptype"].append(sb.get("sp_type", ""))
        cols["exofop_tois"].append("; ".join(ex.get("tois", [])))
        cols["exofop_ctois"].append("; ".join(ex.get("ctois", [])))
        cols["exofop_followup"].append(
            f"imaging {ex.get('n_imaging', 0)}, spectra {ex.get('n_spectroscopy', 0)}, "
            f"time series {ex.get('n_time_series', 0)}" if ex else "")
        cols["exofop_companions"].append("; ".join(ex.get("companions", [])))
        cols["n_variables_nearby"].append(len(s.get("variables", [])))
        cols["period_matches"].append("; ".join(t for _, t in matches))
        cols["n_blend_suspects"].append(len(suspects))
        cols["min_blend_amp_pct"].append(f"{suspects[0][0] * 100:.1f}" if suspects else "")
        cols["blend_suspects"].append("; ".join(t for _, t in suspects[:MAX_SUSPECTS]))
        cols["enrich_flags"].append(";".join(flags))
    for k, v in cols.items():
        df[k] = v

    out = Path(args.out) if args.out else src.parent / "enriched.csv"
    df.to_csv(out, index=False)
    logging.info(f"wrote {len(df)} rows to {out}")
    for _, r in df[df["enrich_flags"] != ""].iterrows():
        logging.info(f"  {r['target']} #{r['candidate']} P={r['period_d']}: {r['enrich_flags']}"
                     + (f" | {r['period_matches']}" if r["period_matches"] else ""))


if __name__ == "__main__":
    main()
