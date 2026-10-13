"""
tess_target_finder.py

Builds a candidate list of dwarf stars for a blind TESS transit search.

Steps:
  1. Query the TESS Input Catalog (TIC) for small, reasonably bright dwarfs
     with low flux contamination, in one patch of sky.
  2. Count TESS sectors per star with tess-point and keep low-coverage stars.
  3. Remove stars that already have a TOI, a confirmed planet, or are known
     TESS eclipsing binaries.
  4. Add Gaia DR3 RUWE (binarity indicator) and drop likely binaries.
  5. Flag whether each star has SPOC 2-min/20-s light curves (i.e. whether the
     official 2-min pipeline will have searched it).
  6. Save a ranked CSV, plus a plain target list ("TIC <number>" per line) that
     run_parallel.sh and find_planets.py read directly.

Requirements:
    pip install astroquery pyvo tess-point pandas numpy
"""

import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astroquery.mast import Catalogs, Observations
from astroquery.vizier import Vizier
from tess_stars2px import tess_stars2px_function_entry

# ---------------------------------------------------------------- settings --
FIELD_NAME = "Aries"     # label recorded in the target file
FIELD_RA = 40.0          # field centre, degrees. Near the ecliptic (few sectors), 40 deg
FIELD_DEC = 15.0         # from the Galactic plane (little crowding), between K2 campaigns
FIELD_RADIUS = 8.5       # degrees; larger = slower queries

TMAG_RANGE = [8, 11]           # TESS magnitude. Injection-recovery on the first field
                               # (T ~ 11-13) gave 50% completeness only at ~3-4 R_earth, so
                               # small planets need brighter stars. Bright dwarfs are rare:
                               # expect ~10x fewer stars than T 8-13, so widen FIELD_RADIUS
RADIUS_RANGE = [0.1, 0.8]      # solar radii
TEFF_RANGE = [3500, 5500]      # kelvin
LOGG_MIN = 4.0
CONTRATIO_MAX = 0.2
RUWE_MAX = 1.4
MIN_SECTORS = 1
MAX_SECTORS = 6

# tess-point also lists sectors scheduled in future. Only count sectors that
# should already have happened (approx: sector 1 began 2018-07-25, ~27.4 d each).
LATEST_SECTOR = int((pd.Timestamp.today() - pd.Timestamp("2018-07-25")).days / 27.4)

HERE = Path(__file__).resolve().parent
OUTPUT_CSV = "targets/aries.csv"        # relative to this folder
OUTPUT_TXT = "targets/aries.txt"        # target list for pipeline.sh / run_parallel.sh

TOI_URL = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
NEA_URL = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
           "query=select+distinct+tic_id+from+pscomppars+where+tic_id+is+not+null"
           "&format=csv")
TESS_EB_VIZIER = "J/ApJS/258/16"   # Prsa et al. (2022) TESS eclipsing binaries

failed_checks = []


GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap"   # ESA Gaia archive


def gaia_query(adql):
    """Run an ADQL query on the Gaia archive (synchronous; up to 2000 rows) and
    return an astropy Table. pyvo, not astroquery.gaia, which contacts the
    archive on import to print its news."""
    import pyvo
    return pyvo.dal.TAPService(GAIA_TAP).run_sync(adql).to_table()


def gaia_text(col):
    """Gaia source ids as text (None where missing): a 19-digit id loses its last
    digits as a float, which silently breaks RUWE lookups."""
    def one(v):
        if v is None or (isinstance(v, float) and np.isnan(v)) or np.ma.is_masked(v):
            return None
        v = str(v).strip().split(".")[0]
        return v if v.isdigit() else None
    return col.map(one)


def safe(func, name, default):
    """Run a lookup; if a service is down, warn loudly instead of crashing."""
    try:
        return func()
    except Exception as exc:  # noqa: BLE001
        print(f"  WARNING: {name} failed ({exc}). This check was SKIPPED.")
        failed_checks.append(name)
        return default


# ------------------------------------------------------------------ steps --
def query_tic():
    print("Querying TIC...")
    tic = Catalogs.query_criteria(
        coordinates=f"{FIELD_RA} {FIELD_DEC}",
        radius=FIELD_RADIUS * u.deg,
        catalog="TIC",
        Tmag=TMAG_RANGE,
        rad=RADIUS_RANGE,
        Teff=TEFF_RANGE,
        contratio=[0, CONTRATIO_MAX],
        lumclass="DWARF",
    )
    df = tic.to_pandas()
    df = df[df["logg"] >= LOGG_MIN].copy()
    df["ID"] = pd.to_numeric(df["ID"]).astype("int64")
    df["GAIA"] = gaia_text(df["GAIA"])            # 19-digit ids don't survive a float
    cols = ["ID", "ra", "dec", "Tmag", "Teff", "rad", "mass", "logg",
            "contratio", "GAIA"]
    print(f"  {len(df)} stars pass the stellar cuts")
    return df[cols].reset_index(drop=True)


def add_sector_counts(df):
    print("Counting TESS sectors (tess-point)...")
    out_id, _, _, out_sec, *_ = tess_stars2px_function_entry(
        df["ID"].values, df["ra"].values, df["dec"].values)
    sec = pd.DataFrame({"ID": out_id, "sector": out_sec})
    # Drop non-observations (-1) and occasional bogus values from tess-point.
    sec = sec[(sec["sector"] > 0) & (sec["sector"] <= LATEST_SECTOR)]
    grouped = sec.groupby("ID")["sector"]
    df["n_sectors"] = df["ID"].map(grouped.nunique()).fillna(0).astype(int)
    df["sectors"] = df["ID"].map(
        grouped.apply(lambda s: ",".join(str(x) for x in sorted(set(s)))))
    print(f"  sectors observed so far (up to ~{LATEST_SECTOR}) per star:")
    print("  " + df["n_sectors"].value_counts().sort_index()
          .to_string().replace("\n", "\n  "))
    df = df[df["n_sectors"].between(MIN_SECTORS, MAX_SECTORS)].copy()
    print(f"  {len(df)} stars have {MIN_SECTORS}-{MAX_SECTORS} sectors")
    return df


def toi_tics():
    toi = pd.read_csv(TOI_URL, comment="#")
    return set(pd.to_numeric(toi["TIC ID"], errors="coerce").dropna().astype("int64"))


def confirmed_tics():
    nea = pd.read_csv(NEA_URL)
    ids = nea["tic_id"].astype(str).str.replace("TIC", "").str.strip()
    return set(pd.to_numeric(ids, errors="coerce").dropna().astype("int64"))


def eb_tics():
    Vizier.ROW_LIMIT = -1
    table = Vizier.get_catalogs(TESS_EB_VIZIER)[0].to_pandas()
    return set(pd.to_numeric(table["TIC"], errors="coerce").dropna().astype("int64"))


def remove_known(df):
    print("Removing known TOIs, planets and eclipsing binaries...")
    known = set()
    known |= safe(toi_tics, "TOI list (ExoFOP)", set())
    known |= safe(confirmed_tics, "Confirmed planets (NASA Exoplanet Archive)", set())
    known |= safe(eb_tics, "TESS eclipsing binary catalogue (VizieR)", set())
    print(f"  {len(known)} known IDs loaded")
    before = len(df)
    df = df[~df["ID"].isin(known)].copy()
    print(f"  removed {before - len(df)}, {len(df)} remain")
    return df


def add_ruwe(df):
    print("Fetching Gaia DR3 RUWE...")
    ids = df["GAIA"].dropna().unique()

    def fetch():
        frames = []
        for chunk in np.array_split(ids, len(ids) // 500 + 1):
            if len(chunk) == 0:
                continue
            q = ("SELECT source_id, ruwe FROM gaiadr3.gaia_source "
                 f"WHERE source_id IN ({','.join(map(str, chunk))})")
            r = gaia_query(q).to_pandas()
            r.columns = [c.lower() for c in r.columns]
            r["source_id"] = r["source_id"].astype("int64").astype(str)
            frames.append(r)
        return pd.concat(frames) if frames else pd.DataFrame(columns=["source_id", "ruwe"])

    r = safe(fetch, "Gaia RUWE", pd.DataFrame(columns=["source_id", "ruwe"]))
    df["ruwe"] = df["GAIA"].map(dict(zip(r["source_id"], r["ruwe"])))
    # Keep stars with no RUWE (TIC uses Gaia DR2 IDs; a few don't match DR3).
    df = df[(df["ruwe"] <= RUWE_MAX) | df["ruwe"].isna()].copy()
    print(f"  {len(df)} remain after RUWE cut ({df['ruwe'].isna().sum()} without RUWE)")
    return df


def flag_short_cadence(df):
    print("Checking for SPOC 2-min/20-s light curves...")
    ids = df["ID"].astype(str).tolist()

    def fetch():
        found = set()
        for i in range(0, len(ids), 200):
            obs = Observations.query_criteria(
                obs_collection="TESS", dataproduct_type="timeseries",
                target_name=ids[i:i + 200])
            found.update(str(t) for t in obs["target_name"])
        return found

    found = safe(fetch, "SPOC short-cadence check", set())
    df["has_spoc_2min"] = df["ID"].astype(str).isin(found)
    return df


def add_duration_hint(df):
    # Max (central) transit duration at P = 10 d, from the rule of thumb.
    rho = df["mass"] / df["rad"] ** 3          # solar units
    df["tmax_10d_hours"] = (13 * (10 / 365) ** (1 / 3) * rho ** (-1 / 3)).round(2)
    return df


def criteria_lines(counts=None):
    """The selection criteria (and star counts per step, if given), as '#' lines
    for the head of a target file, so the list records how it was made."""
    lines = [
        f"# Selection (tess_target_finder.py, run {pd.Timestamp.today():%Y-%m-%d}):",
        f"#   field        {FIELD_NAME}: RA {FIELD_RA} deg, Dec {FIELD_DEC} deg, radius {FIELD_RADIUS} deg",
        f"#   TESS mag     {TMAG_RANGE[0]} to {TMAG_RANGE[1]}",
        f"#   star         dwarf (TIC lumclass), radius {RADIUS_RANGE[0]}-{RADIUS_RANGE[1]} R_sun, "
        f"Teff {TEFF_RANGE[0]}-{TEFF_RANGE[1]} K, log g >= {LOGG_MIN}",
        f"#   blending     TIC contamination ratio <= {CONTRATIO_MAX}; Gaia DR3 RUWE <= {RUWE_MAX} "
        f"(stars without RUWE kept)",
        f"#   coverage     {MIN_SECTORS}-{MAX_SECTORS} TESS sectors observed (tess-point, sectors <= "
        f"{LATEST_SECTOR})",
        "#   removed      known TOIs (ExoFOP), confirmed planets (NASA Exoplanet Archive),",
        "#                TESS eclipsing binaries (Prsa et al. 2022)",
    ]
    if counts:
        lines.append("#   stars left   " + " -> ".join(f"{n} {label}" for label, n in counts))
    if failed_checks:
        lines.append("#   SKIPPED      " + "; ".join(failed_checks) + " (service failed: list less clean)")
    return lines


def write_target_list(df, path, counts=None):
    """One "TIC <number>" per line, least blended first, so the most trustworthy
    stars are searched first if a run is stopped early. '#' starts a comment.
    The header records the selection criteria."""
    lines = [
        f"# {FIELD_NAME} field: {len(df)} targets from {OUTPUT_CSV}.",
        *criteria_lines(counts),
        "# Order: contamination ratio (least blended first), so the most trustworthy",
        "# stars are searched first if the run is stopped early.",
        "# Comment after each name: Tmag, Teff, radius (R_sun), contamination ratio, sectors,",
        "# and '2min' if SPOC 2-min light curves exist (already searched by the official pipeline).",
        f"# Run with:  ./pipeline.sh {path} results/<name> 8",
        "",
    ]
    for r in df.sort_values("contratio").itertuples():
        lines.append(f"{r.target:<15} # T={r.Tmag:.1f}  Teff={r.Teff:.0f}  R={r.rad:.2f}  "
                     f"contam={r.contratio:.3f}  sectors={r.sectors}"
                     + ("  2min" if r.has_spoc_2min else ""))
    (HERE / path).write_text("\n".join(lines) + "\n")


# ---------------------------------------------------- whole-sky catalogue --
# The best targets are a property of each star, not of a patch of sky: the
# catalogue holds every star in the sky that passes the stellar cuts, scored by
# how well a small planet could be found around it, and batches hand out the
# best stars not yet searched. Stored in tess.db (catalog, catalog_tiles, batches).
TILE_DEG = 15                    # sky queried in RA/Dec boxes of this size...
TILE_WORKERS = 4                 # ...this many at a time (MAST allows a few in parallel)
TILE_TRIES = 3
SCORE_RP = 2.0                   # score = expected SNR of this planet (Earth radii)...
SCORE_P = 10.0                   # ...on this period (days)
USABLE_DAYS = 24.0               # usable days of data per sector


def noise_1h_ppm(tmag):
    """Approximate TESS noise in 1 hour (ppm): photon/background noise scaling with
    brightness plus a ~60 ppm systematic floor. A rough model, for ranking only."""
    return np.sqrt((230 * 10 ** (0.2 * (np.asarray(tmag, float) - 10))) ** 2 + 60 ** 2)


def expected_snr(df, rp=SCORE_RP, period=SCORE_P):
    """Expected SNR of a central transit of an rp planet on `period` days, from the
    star's T mag, radius, mass and number of sectors. Ranks stars; injection tests
    put 50% completeness near an expected SNR of 12."""
    R = df["rad"].astype(float)
    M = df["mass"].astype(float).fillna(R)                      # dwarfs: M ~ R in solar units
    k = rp / 109.1 / R
    a = 215.032 * (M * (period / 365.25) ** 2) ** (1 / 3) / R    # a / R*
    t14_h = period / np.pi * np.arcsin(np.clip((1 + k) / a, 0, 1)) * 24
    n_tr = df["n_sectors"].astype(float) * USABLE_DAYS / period
    return k ** 2 * 1e6 / noise_1h_ppm(df["tmag"]) * np.sqrt(t14_h * n_tr)


def db():
    """tess.db, with its schema (catalogue tables included)."""
    import tessdb
    return tessdb.connect()


def query_tile(ra0, dec0):
    """One RA/Dec box of the TIC with the stellar cuts, as a DataFrame."""
    import time
    for attempt in range(1, TILE_TRIES + 1):
        try:
            tab = Catalogs.query_criteria(
                catalog="TIC", ra=[ra0, ra0 + TILE_DEG], dec=[dec0, dec0 + TILE_DEG],
                Tmag=TMAG_RANGE, rad=RADIUS_RANGE, Teff=TEFF_RANGE,
                contratio=[0, CONTRATIO_MAX], lumclass="DWARF")
            df = tab.to_pandas()
            df["GAIA"] = gaia_text(df["GAIA"])
            return df[df["logg"] >= LOGG_MIN]
        except Exception as e:
            if attempt == TILE_TRIES:
                raise
            print(f"  tile RA {ra0} Dec {dec0}: attempt {attempt} failed ({type(e).__name__}); retrying")
            time.sleep(30 * attempt)


def build_catalog(tile_limit=None):
    """Query the whole sky tile by tile (resumable), then add sectors, known-system
    flags, Gaia RUWE and the score to every star in the catalog table."""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from astropy.coordinates import SkyCoord
    con = db()
    done = {r[0] for r in con.execute("SELECT tile FROM catalog_tiles")}
    tiles = [(ra, dec) for dec in range(-90, 90, TILE_DEG) for ra in range(0, 360, TILE_DEG)
             if f"{ra}:{dec}" not in done]
    if tile_limit:
        tiles = tiles[:tile_limit]                       # testing only
    print(f"Whole-sky catalogue: {len(done)} tiles already done, {len(tiles)} to query "
          f"({TILE_WORKERS} at a time)")
    cols = ["ID", "ra", "dec", "Tmag", "Teff", "rad", "mass", "logg", "contratio", "GAIA"]
    now = pd.Timestamp.now().isoformat(timespec="seconds")
    with ThreadPoolExecutor(TILE_WORKERS) as pool:
        futures = {pool.submit(query_tile, ra, dec): (ra, dec) for ra, dec in tiles}
        for i, f in enumerate(as_completed(futures), 1):
            ra, dec = futures[f]
            try:
                df = f.result()
            except Exception as e:
                print(f"  tile RA {ra} Dec {dec} FAILED ({type(e).__name__}); rerun to retry it")
                continue
            with con:
                for r in df[cols].itertuples(index=False):
                    con.execute("INSERT INTO catalog (tic, ra, dec, tmag, teff, rad, mass, logg, contratio,"
                                " gaia, updated) VALUES (?,?,?,?,?,?,?,?,?,?,?)"
                                " ON CONFLICT (tic) DO UPDATE SET tmag=excluded.tmag, updated=excluded.updated",
                                (str(int(r.ID)), r.ra, r.dec, r.Tmag, r.Teff, r.rad,
                                 None if pd.isna(r.mass) else r.mass, r.logg, r.contratio, r.GAIA, now))
                con.execute("INSERT OR REPLACE INTO catalog_tiles VALUES (?,?,?)", (f"{ra}:{dec}", len(df), now))
            if i % 12 == 0 or i == len(tiles):
                print(f"  {i}/{len(tiles)} tiles; catalogue now "
                      f"{con.execute('SELECT count(*) FROM catalog').fetchone()[0]} stars")
    n_tiles = len(range(-90, 90, TILE_DEG)) * len(range(0, 360, TILE_DEG))
    left = n_tiles - con.execute("SELECT count(*) FROM catalog_tiles").fetchone()[0]
    cat = pd.read_sql("SELECT * FROM catalog", con)
    print(f"{len(cat)} stars pass the stellar cuts across the sky")
    # sectors observed so far (tess-point runs locally)
    out_id, _, _, out_sec, *_ = tess_stars2px_function_entry(cat["tic"].astype("int64").values,
                                                             cat["ra"].values, cat["dec"].values)
    sec = pd.DataFrame({"tic": out_id.astype(str), "sector": out_sec})
    sec = sec[(sec["sector"] > 0) & (sec["sector"] <= LATEST_SECTOR)]
    g = sec.groupby("tic")["sector"]
    cat["n_sectors"] = cat["tic"].map(g.nunique()).fillna(0).astype(int)
    cat["sectors"] = cat["tic"].map(g.apply(lambda x: ",".join(str(v) for v in sorted(set(x)))))
    # coordinates for later choices (crowding, ecliptic)
    c = SkyCoord(cat["ra"].values * u.deg, cat["dec"].values * u.deg)
    cat["glat"] = c.galactic.b.deg
    cat["elat"] = c.barycentrictrueecliptic.lat.deg
    # known systems
    print("Flagging known TOIs, planets and eclipsing binaries...")
    known = {}
    for label, fn in (("TOI", toi_tics), ("planet", confirmed_tics), ("EB", eb_tics)):
        for t in safe(fn, label, set()):
            known.setdefault(str(t), label)
    cat["known"] = cat["tic"].map(known)
    # Gaia RUWE, only where missing (resumable)
    need = cat[cat["ruwe"].isna() & cat["gaia"].notna()]
    ruwe = {}
    if len(need):
        print(f"Fetching Gaia DR3 RUWE for {len(need)} stars...")
        ids = need["gaia"].unique()
        for chunk in [ids[i:i + 500] for i in range(0, len(ids), 500)]:
            q = f"SELECT source_id, ruwe FROM gaiadr3.gaia_source WHERE source_id IN ({','.join(chunk)})"
            res = safe(lambda: gaia_query(q).to_pandas(), "Gaia RUWE", pd.DataFrame())
            if len(res):
                res.columns = [x.lower() for x in res.columns]
                ruwe.update(zip(res["source_id"].astype("int64").astype(str), res["ruwe"]))
        cat.loc[need.index, "ruwe"] = need["gaia"].map(ruwe)
    cat["score"] = expected_snr(cat)
    with con:
        for r in cat.itertuples(index=False):
            con.execute("UPDATE catalog SET n_sectors=?, sectors=?, glat=?, elat=?, known=?, ruwe=?, score=?"
                        " WHERE tic=?", (int(r.n_sectors), r.sectors, r.glat, r.elat,
                                         None if pd.isna(r.known) else r.known,
                                         None if pd.isna(r.ruwe) else float(r.ruwe), float(r.score), r.tic))
    eligible = cat[(cat["n_sectors"].between(MIN_SECTORS, MAX_SECTORS)) & cat["known"].isna()
                   & ((cat["ruwe"] <= RUWE_MAX) | cat["ruwe"].isna())]
    print(f"\nCatalogue: {len(cat)} stars; {len(eligible)} eligible ({MIN_SECTORS}-{MAX_SECTORS} sectors, "
          f"not known, RUWE <= {RUWE_MAX}); median score {eligible['score'].median():.1f}; "
          f"{(eligible['score'] >= 12).sum()} with expected SNR >= 12 for a {SCORE_RP} R_earth, "
          f"{SCORE_P:.0f} d planet")
    if left > 0:
        print(f"WARNING: {left} tile(s) failed; rerun --build-catalog to fill them")


def make_batch(n, name):
    """Write the next batch: the n highest-scoring eligible stars that haven't been
    searched (any run in tess.db) or handed out in an earlier batch."""
    con = db()
    if con.execute("SELECT 1 FROM batches WHERE name = ?", (name,)).fetchone():
        raise SystemExit(f"batch {name} already exists")
    cat = pd.read_sql("SELECT * FROM catalog WHERE score IS NOT NULL", con)
    if cat.empty:
        raise SystemExit("no catalogue yet: run with --build-catalog first")
    searched = {r[0] for r in con.execute("SELECT DISTINCT star FROM searches")}
    pool = cat[cat["n_sectors"].between(MIN_SECTORS, MAX_SECTORS) & cat["known"].isna()
               & ((cat["ruwe"] <= RUWE_MAX) | cat["ruwe"].isna()) & cat["batch"].isna()
               & ~("TIC " + cat["tic"]).isin(searched)]
    df = pool.sort_values("score", ascending=False).head(n).copy()
    if df.empty:
        raise SystemExit("no eligible stars left")
    df = df.rename(columns={"tic": "ID", "tmag": "Tmag", "teff": "Teff", "gaia": "GAIA"})
    df = flag_short_cadence(df)
    df.insert(0, "target", "TIC " + df["ID"].astype(str))
    csv_path, txt_path = f"targets/{name}.csv", f"targets/{name}.txt"
    df.to_csv(HERE / csv_path, index=False)
    built = con.execute("SELECT max(done) FROM catalog_tiles").fetchone()[0]
    lines = [
        f"# Batch {name}: the {len(df)} highest-scoring unsearched stars in the whole-sky catalogue",
        f"# (built {built}), written {pd.Timestamp.today():%Y-%m-%d} by tess_target_finder.py --batch.",
        f"#   TESS mag     {TMAG_RANGE[0]} to {TMAG_RANGE[1]}",
        f"#   star         dwarf (TIC lumclass), radius {RADIUS_RANGE[0]}-{RADIUS_RANGE[1]} R_sun, "
        f"Teff {TEFF_RANGE[0]}-{TEFF_RANGE[1]} K, log g >= {LOGG_MIN}",
        f"#   blending     TIC contamination ratio <= {CONTRATIO_MAX}; Gaia DR3 RUWE <= {RUWE_MAX}",
        f"#   coverage     {MIN_SECTORS}-{MAX_SECTORS} sectors (tess-point, sectors <= {LATEST_SECTOR})",
        "#   removed      known TOIs, confirmed planets, TESS eclipsing binaries; stars already",
        "#                searched (tess.db) or in an earlier batch",
        f"#   ranking      expected SNR of a {SCORE_RP} R_earth planet on a {SCORE_P:.0f} d orbit "
        f"(score {df['score'].min():.1f}-{df['score'].max():.1f})",
        f"#   eligible     {len(pool)} stars before this batch",
        "# Order: highest score first. Comment: score, Tmag, Teff, radius, contamination, sectors, '2min'.",
        f"# Prefetch (optional, any time MAST is healthy):  .venv/bin/python prefetch.py {txt_path} results/{name}",
        f"# Run with:  ./pipeline.sh {txt_path} results/{name} 8",
        "",
    ]
    for r in df.itertuples():
        lines.append(f"{r.target:<15} # score={r.score:.1f}  T={r.Tmag:.1f}  Teff={r.Teff:.0f}  R={r.rad:.2f}  "
                     f"contam={r.contratio:.3f}  sectors={r.sectors}" + ("  2min" if r.has_spoc_2min else ""))
    (HERE / txt_path).write_text("\n".join(lines) + "\n")
    with con:
        con.executemany("UPDATE catalog SET batch = ? WHERE tic = ?", [(name, str(t)) for t in df["ID"]])
        con.execute("INSERT INTO batches VALUES (?,?,?,?,?)",
                    (name, pd.Timestamp.now().isoformat(timespec="seconds"), len(df),
                     "\n".join(lines[2:11]), txt_path))
    print(f"Batch {name}: {len(df)} stars (score {df['score'].max():.1f} to {df['score'].min():.1f}; "
          f"{int(df['has_spoc_2min'].sum())} with SPOC 2-min) -> {txt_path}")


def release_batch(name):
    """Return an unused batch's stars to the pool (e.g. after changing the selection
    criteria) and delete its target files. Refused once the batch has been run or
    started: its results/<name> folder holds search output (a prefetch alone doesn't
    count), or any of its stars was searched."""
    con = db()
    row = con.execute("SELECT targets_file FROM batches WHERE name = ?", (name,)).fetchone()
    if not row:
        raise SystemExit(f"no batch named {name}")
    out = HERE / "results" / name
    started = [f.name for f in out.glob("*") if f.name != "prefetch.jsonl"] if out.is_dir() else []
    if started:
        raise SystemExit(f"batch {name} has been run or started (results/{name} has "
                         f"{', '.join(sorted(started)[:3])}): not released")
    n_searched = con.execute("SELECT count(DISTINCT s.star) FROM searches s JOIN catalog c "
                             "ON s.star = 'TIC ' || c.tic WHERE c.batch = ?", (name,)).fetchone()[0]
    if n_searched:
        raise SystemExit(f"batch {name}: {n_searched} of its stars have been searched: not released")
    with con:
        n = con.execute("UPDATE catalog SET batch = NULL WHERE batch = ?", (name,)).rowcount
        con.execute("DELETE FROM batches WHERE name = ?", (name,))
    removed = []
    for f in (row[0], str(Path(row[0]).with_suffix(".csv"))):
        if f and (HERE / f).exists():
            (HERE / f).unlink()
            removed.append(f)
    if (out / "prefetch.jsonl").exists():          # a prefetch only: its manifest goes too
        (out / "prefetch.jsonl").unlink()
        out.rmdir()
        removed.append(f"results/{name}/prefetch.jsonl (its downloads stay in the lightkurve cache)")
    print(f"Released batch {name}: {n} stars back in the pool"
          + (f"; deleted {', '.join(removed)} (commit the removal if they were in git)" if removed else ""))


# ------------------------------------------------------------------- main --
def main():
    import argparse
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser(description="Choose TESS targets: a cone (default), or the whole-sky "
                                             "catalogue (--build-catalog, then --batch N --name NAME; --release NAME undoes an unused batch).")
    ap.add_argument("--build-catalog", action="store_true", help="query the whole sky into tess.db (resumable)")
    ap.add_argument("--batch", type=int, metavar="N", help="write the next N unsearched stars from the catalogue")
    ap.add_argument("--name", help="batch name (for --batch)")
    ap.add_argument("--release", metavar="NAME", help="return an unused batch's stars to the pool")
    ap.add_argument("--tile-limit", type=int, help=argparse.SUPPRESS)          # testing
    args = ap.parse_args()
    if args.build_catalog:
        build_catalog(args.tile_limit)
        return
    if args.release:
        release_batch(args.release)
        return
    if args.batch:
        if not args.name:
            ap.error("--batch needs --name")
        make_batch(args.batch, args.name)
        return
    df = query_tic()
    counts = [("pass the stellar cuts", len(df))]
    labels = {"add_sector_counts": "with the sector range", "remove_known": "after removing known",
              "add_ruwe": "after the RUWE cut"}
    for step in (add_sector_counts, remove_known, add_ruwe,
                 flag_short_cadence, add_duration_hint):
        if df.empty:
            print("\nNo stars left - loosen the settings or try another field.")
            return
        df = step(df)
        if step.__name__ in labels:
            counts.append((labels[step.__name__], len(df)))

    # Rank: FFI-only stars first (never searched by the 2-min pipeline),
    # then smallest stars (deepest transits), then brightest.
    df = df.sort_values(["has_spoc_2min", "rad", "Tmag"]).reset_index(drop=True)
    df.insert(0, "target", "TIC " + df["ID"].astype(str))   # name lightkurve/MAST accept
    df.to_csv(HERE / OUTPUT_CSV, index=False)
    write_target_list(df, OUTPUT_TXT, counts)
    print(f"\nSaved {len(df)} targets to {OUTPUT_CSV} and {OUTPUT_TXT}")
    if failed_checks:
        print("WARNING - these checks were skipped, so treat the list with caution:")
        for name in failed_checks:
            print(f"  - {name}")


if __name__ == "__main__":
    main()
