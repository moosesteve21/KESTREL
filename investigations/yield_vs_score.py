#!/usr/bin/env python3
"""Yield of the whole-sky batches against the ranking score and other star properties.

For every star searched in the b... batches: its score (expected SNR of a 2 R_earth, 10 d
planet), Tmag, Teff, sector count and whether it has SPOC 2-min data (targets/b*.csv), and
its outcome: any signal, a clean signal, a kept candidate (on the watchlist), or a clean
signal that was already known (TOI, CTOI or SPOC TCE). Prints rates per bin with Poisson
errors (small numbers: read the error bars).

Usage: .venv/bin/python investigations/yield_vs_score.py [--db db/tess.db]
"""
import argparse
import glob
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load(db):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    cat = pd.read_sql("SELECT tic, tmag, teff, rad, n_sectors, score, batch FROM catalog "
                      "WHERE batch LIKE 'b%'", con)
    cat["star"] = "TIC " + cat["tic"].astype(str)
    runs = pd.read_sql("SELECT r.name AS run, s.star FROM searches s JOIN runs r USING (run_id) "
                       "WHERE r.name LIKE 'b%'", con)
    cand = pd.read_sql("SELECT r.name AS run, c.star, c.status, c.enrich_flags FROM candidates c "
                       "JOIN runs r USING (run_id) WHERE r.name LIKE 'b%' AND c.status != 'same_period'", con)
    watch = {r[0] for r in con.execute("SELECT star FROM watchlist")}
    two_min = {}
    for f in glob.glob(str(ROOT / "targets" / "b*.csv")):
        t = pd.read_csv(f, usecols=lambda c: c in ("target", "has_spoc_2min"))
        if "has_spoc_2min" in t:
            two_min.update(dict(zip(t["target"], t["has_spoc_2min"].astype(bool))))
    df = cat[cat["star"].isin(set(runs["star"]))].copy()
    df["two_min"] = df["star"].map(two_min)
    df["signal"] = df["star"].isin(set(cand["star"]))
    clean = cand[cand["status"] == "clean"]
    df["clean"] = df["star"].isin(set(clean["star"]))
    known = clean[clean["enrich_flags"].fillna("").str.contains("has_toi|has_ctoi|has_spoc_tce")]
    df["known"] = df["star"].isin(set(known["star"]))
    df["kept"] = df["star"].isin(watch)
    return df


def table(df, col, bins, labels):
    df = df.assign(bin=pd.cut(df[col], bins=bins, labels=labels, include_lowest=True))
    print(f"\n{'by ' + col:<18} {'stars':>6} {'signal %':>9} {'clean %':>9} {'known':>6} {'kept':>5} "
          f"{'kept per 1000':>15}")
    for b, g in df.groupby("bin", observed=True):
        n, k = len(g), int(g["kept"].sum())
        print(f"{str(b):<18} {n:>6} {100 * g['signal'].mean():>8.1f}% {100 * g['clean'].mean():>8.1f}% "
              f"{int(g['known'].sum()):>6} {k:>5} {1000 * k / n:>8.1f} +- {1000 * np.sqrt(max(k, 1)) / n:.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=str(ROOT / "db" / "tess.db"))
    df = load(ap.parse_args().db)
    print(f"{len(df)} searched stars from the whole-sky batches; {int(df['signal'].sum())} with a signal, "
          f"{int(df['clean'].sum())} clean, {int(df['known'].sum())} clean but already known, "
          f"{int(df['kept'].sum())} kept candidates")
    q = np.unique(np.quantile(df["score"], [0, 0.2, 0.4, 0.6, 0.8, 1]))
    table(df, "score", q, [f"{a:.1f}-{b:.1f}" for a, b in zip(q[:-1], q[1:])])
    if df["two_min"].notna().any():
        d2 = df[df["two_min"].notna()].assign(two=lambda x: x["two_min"].astype(int))
        table(d2, "two", [-0.5, 0.5, 1.5], ["FFI only", "SPOC 2-min"])
    table(df, "n_sectors", [0.5, 1.5, 2.5, 3.5, 4.5, 6.5], ["1", "2", "3", "4", "5-6"])
    table(df, "teff", [2000, 3900, 4600, 5200, 6000], ["<3900 K (M)", "3900-4600", "4600-5200", "5200-6000"])
    table(df, "tmag", [0, 9, 10, 10.5, 11, 14], ["<9", "9-10", "10-10.5", "10.5-11", ">11"])


if __name__ == "__main__":
    main()
