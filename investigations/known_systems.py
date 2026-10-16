#!/usr/bin/env python3
"""known_systems.py - the known-system recovery sample for the paper.

The whole-sky batches left out every catalogue star that was a known system when the
catalogue was built (2026-10-11): a TOI host, a confirmed-planet host or a TESS eclipsing
binary (catalog.known; the first label found wins, in that order). This picks exactly
those stars that pass every other batch rule (1-6 sectors, RUWE <= 1.4 or unknown, a
score; the stellar cuts were applied when the catalogue was built), so the frozen
pipeline can be run on them blind and scored against what is known:
  - recovery: confirmed and known planets (TOI disposition CP/KP, or the Exoplanet Archive);
  - rejection: TOIs dispositioned FP/FA and the TESS eclipsing binaries;
  - PC TOIs (unconfirmed) are reported separately.

It saves the "truth" as it stands today, so the sample can be reproduced exactly:
  investigations/known_systems/toi_<date>.csv.gz        ExoFOP TOI list (full)
  investigations/known_systems/nea_<date>.csv.gz        Exoplanet Archive pscomppars (sample stars)
  investigations/known_systems/eb_<date>.csv.gz         Prsa et al. 2022 TESS EBs (sample stars)
  investigations/known_systems/truth.csv                one row per known object on a sample star
  targets/<name>.txt                                    the target list (score order)

Usage: .venv/bin/python investigations/known_systems.py [--name k001]
Run the list later with: ./pipeline.sh targets/k001.txt results/k001 8 "known-system recovery k001"
"""
import argparse
import gzip
import io
import sqlite3
import subprocess
import sys
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import tess_target_finder as ttf  # noqa: E402  (the batch rules and source URLs, unchanged)

OUT = ROOT / "investigations" / "known_systems"
NEA = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,tic_id,pl_orbper,"
       "pl_tranmid,pl_trandur,pl_trandep,pl_rade,tran_flag,disc_facility+from+pscomppars"
       "+where+tic_id+is+not+null&format=csv")


def fetch(url):
    with urllib.request.urlopen(url, timeout=300) as r:
        return r.read()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="k001")
    name = ap.parse_args().name
    today = date.today().isoformat()
    OUT.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(f"file:{ROOT / 'db' / 'tess.db'}?mode=ro", uri=True)
    cat = pd.read_sql("SELECT * FROM catalog WHERE score IS NOT NULL", con)
    s = cat[cat["n_sectors"].between(ttf.MIN_SECTORS, ttf.MAX_SECTORS) & cat["known"].notna()
            & ((cat["ruwe"] <= ttf.RUWE_MAX) | cat["ruwe"].isna())].sort_values("score", ascending=False)
    tics = set(s["tic"].astype("int64"))
    print(f"sample: {len(s)} stars ({s['known'].value_counts().to_dict()})")

    toi_raw = fetch(ttf.TOI_URL)
    (OUT / f"toi_{today}.csv.gz").write_bytes(gzip.compress(toi_raw))
    toi = pd.read_csv(io.BytesIO(toi_raw), comment="#")
    toi = toi[pd.to_numeric(toi["TIC ID"], errors="coerce").isin(tics)]
    nea = pd.read_csv(io.BytesIO(fetch(NEA)))
    nea["tic"] = pd.to_numeric(nea["tic_id"].astype(str).str.replace("TIC", "").str.strip(), errors="coerce")
    nea = nea[nea["tic"].isin(tics)]
    nea.to_csv(OUT / f"nea_{today}.csv.gz", index=False)
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    eb = Vizier.get_catalogs(ttf.TESS_EB_VIZIER)[0].to_pandas()
    eb = eb[pd.to_numeric(eb["TIC"], errors="coerce").isin(tics)]
    eb.to_csv(OUT / f"eb_{today}.csv.gz", index=False)

    rows = []
    for _, r in toi.iterrows():
        disp = str(r.get("TFOPWG Disposition") or "").strip() or str(r.get("TESS Disposition") or "").strip()
        rows.append(dict(tic=int(r["TIC ID"]), source="TOI", id=f"TOI-{r['TOI']}", disposition=disp,
                         period_d=r.get("Period (days)"), t0_bjd=r.get("Epoch (BJD)"),
                         duration_h=r.get("Duration (hours)"), depth_ppm=r.get("Depth (ppm)"),
                         rp_earth=r.get("Planet Radius (R_Earth)")))
    for _, r in nea.iterrows():
        rows.append(dict(tic=int(r["tic"]), source="NEA", id=r["pl_name"], disposition="CONFIRMED",
                         period_d=r["pl_orbper"], t0_bjd=r["pl_tranmid"], duration_h=r["pl_trandur"],
                         depth_ppm=(r["pl_trandep"] * 1e4 if pd.notna(r["pl_trandep"]) else None),
                         rp_earth=r["pl_rade"], transiting=r["tran_flag"]))
    pcol = next((c for c in eb.columns if c.lower() in ("per", "period", "p")), None)
    tcol = next((c for c in eb.columns if c.lower() in ("bjd0", "t0", "epoch")), None)
    for _, r in eb.iterrows():
        rows.append(dict(tic=int(r["TIC"]), source="TESS-EB", id=f"EB TIC {int(r['TIC'])}", disposition="EB",
                         period_d=r[pcol] if pcol else None, t0_bjd=r[tcol] if tcol else None))
    truth = pd.DataFrame(rows).sort_values(["tic", "source", "period_d"])
    truth.insert(0, "star", "TIC " + truth["tic"].astype(str))
    truth.drop(columns="tic").to_csv(OUT / "truth.csv", index=False)
    missing = tics - set(truth["tic"])
    print(f"truth: {len(truth)} known objects on {truth['tic'].nunique()} stars; "
          f"{len(missing)} sample stars with no entry today (dropped from a list since the build?)")
    print(truth.groupby("source")["disposition"].value_counts().to_string())

    version = subprocess.run(["git", "-C", str(ROOT), "describe", "--tags", "--always", "--dirty"],
                             capture_output=True, text=True).stdout.strip()
    lines = [f"# Known-system recovery sample {name}: {len(s)} catalogue stars that were known systems when the",
             f"# catalogue was built (TOI {int((s['known'] == 'TOI').sum())}, confirmed planet "
             f"{int((s['known'] == 'planet').sum())}, TESS EB {int((s['known'] == 'EB').sum())}) and pass every other",
             f"# batch rule ({ttf.MIN_SECTORS}-{ttf.MAX_SECTORS} sectors, RUWE <= {ttf.RUWE_MAX}). Written {today} by "
             f"investigations/known_systems.py, code {version}.",
             f"# Truth: investigations/known_systems/truth.csv (TOI list, Exoplanet Archive, TESS EBs as of {today}).",
             f"# Prefetch:  .venv/bin/python prefetch.py targets/{name}.txt results/{name}",
             f"# Run with:  ./pipeline.sh targets/{name}.txt results/{name} 8 \"known-system recovery {name}\"", ""]
    for r in s.itertuples():
        lines.append(f"TIC {r.tic}   # {r.known}, score {r.score:.1f}, T={r.tmag:.2f}, Teff {r.teff:.0f}, "
                     f"{r.n_sectors} sectors")
    (ROOT / "targets" / f"{name}.txt").write_text("\n".join(lines) + "\n")
    print(f"written targets/{name}.txt and {OUT.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
