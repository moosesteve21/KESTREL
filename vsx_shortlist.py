#!/usr/bin/env python3
"""
vsx_shortlist.py - list a run's signals worth a look for the weekly VSX session.

Usage:
    python vsx_shortlist.py results/<name>             # writes results/<name>/vsx_shortlist.csv
    python vsx_shortlist.py results/b0*                # several runs

Reads the run's enriched.csv and enrich_cache/ (no network), and the pixel and
reduction checks' summaries when present, so run it after those checks. A star is listed when its signal looks like a variable star
rather than a transit, and VSX either lacks it or disagrees with TESS:

  signature (any of): sinusoidal_variability or significant_secondary on the
      signal; a same-period row that is a secondary_eclipse, or a phase_curve_residual
      at P <= 3 d (rotation or ellipsoidal; at longer periods it's mostly noise);
      odd_even_mismatch with enough transits (an eclipsing binary found at half its
      period)
  skipped: TOIs and CTOIs, known planet hosts, and signals matching a catalogued
      variable on another star (nearby_variable_period_match: the neighbour's
      variability, already known)
  groups:
      new     no VSX entry at the target (other catalogues that have it are listed:
              Gaia DR3 classes and ZTF entries are often missing from VSX)
      published  not in VSX, or in VSX without this period, but the TESS Ten Thousand
              Catalog (Kostov et al. 2025, 2025ApJS..279...50K) has it as an eclipsing
              binary at this period: submit it as a published star (or a revision)
              with that reference, not as a discovery
      period  in VSX, but without a period, or with one TESS doesn't confirm: no
              match at all, or a match only as a one-day alias (TESS observes
              continuously, so an alias match means the catalogue period is the
              alias - e.g. a catalogued 20.82 d against TESS's 0.95 d)
      (in VSX with the TESS period, or a harmonic of it: not listed)
      blend   the pixel check puts the dip OFF-TARGET and a star there is bright enough
              to cause it (needing < 50% dimming): probably a variable neighbour, named
              in blend_suspects, which may itself be new to VSX. If no star there is
              bright enough, the signal keeps its own group, with a note.
  dropped: signals the reduction check calls reduction_dependent (not in other
      reductions of the same data: most likely systematics, not a variable)

Stars with a vsx, vsx-candidate or submitted note in tess.db are listed with that
note, so the queue isn't repeated. It only makes the list: whether a star goes to
VSX is decided in the weekly session. If the run is already in tess.db, its list is
loaded there too (table vsx_shortlist, view vsx_shortlist_view).
"""
import argparse
import json
import re
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
TARGET_ARCSEC = 5.0             # a catalogue entry this close is the target itself (as enrich)
SIGNATURE_FLAGS = ("sinusoidal_variability", "significant_secondary")
RELATED_FLAGS = ("secondary_eclipse", "phase_curve_residual")
SKIP_ENRICH = ("has_toi", "has_ctoi", "known_planet_period_match", "nearby_variable_period_match")
PHASE_CURVE_MAX_P = 3.0         # days: a phase_curve_residual row counts below this
NOT_EB_FLAGS = ("few_transits", "single_event_dominated")   # odd/even from a handful of dips
BLEND_MAX_DIM = 0.5             # a neighbour needing a deeper eclipse than this can't be the source
HANDLED_TAGS = ("vsx", "vsx-candidate", "submitted")
COLUMNS = ["group", "star", "candidate", "period_d", "depth_ppm", "snr", "signature", "vsx",
           "other_catalogues", "simbad_type", "blend_suspects", "note", "flags", "enrich_flags"]


def flags_of(v):
    return set(filter(None, str(v if isinstance(v, str) else "").split(";")))


def signature(row, related):
    """Why the signal looks like a variable star (empty: it doesn't)."""
    f = flags_of(row["flags"])
    why = [x for x in SIGNATURE_FLAGS if x in f]
    why += sorted({x for r in related for x in flags_of(r["flags"]) if x in RELATED_FLAGS
                   and (x != "phase_curve_residual" or float(row["period_d"]) <= PHASE_CURVE_MAX_P)})
    if "odd_even_mismatch" in f and not f & set(NOT_EB_FLAGS):
        why.append("odd_even_mismatch")
    return why


def handled_notes():
    try:
        import tessdb
        con = sqlite3.connect(tessdb.DB)
        q = f"SELECT star, tag, text FROM notes WHERE tag IN ({','.join('?' * len(HANDLED_TAGS))})"
        out = {}
        for star, tag, text in con.execute(q, HANDLED_TAGS):
            out.setdefault(star, f"{tag}: {text[:80]}")
        return out
    except Exception:
        return {}


def off_target_sources(verdict, depth, tmag):
    """From a pixel verdict's 'OFF-TARGET ...: TIC 1 (T=14.7, 3"); ...' list, the stars
    that could make a dip of this depth: (plausible, a note on the brightest too-faint one).
    A neighbour must dim by depth x 10^(0.4 (T_n - T_target)) or more (more still if it sits
    off the aperture), so above BLEND_MAX_DIM it can't be the source."""
    listed = re.findall(r"TIC (\d+) \(T=([\d.]+), ([\d.]+)\"\)", verdict.split("OFF-TARGET", 1)[1])
    if tmag is None or not np.isfinite(tmag):
        return [f"TIC {t} (T={m}, {d}\")" for t, m, d in listed], ""
    plausible, faint = [], ""
    for tic, m, sep in listed:
        need = depth * 10 ** (0.4 * (float(m) - tmag))
        if need <= BLEND_MAX_DIM:
            plausible.append(f"TIC {tic} (T={m}, {sep}\", needs {need:.1%})")
        elif not faint:
            faint = f"TIC {tic}, T={m}, would need {need:.0%}"
    return plausible, faint


def later_checks(folder):
    """{(star, candidate): (pixel verdict, reduction verdict)} from the run's
    pixel_check*/ and reduction_check*/ summaries, where present."""
    pix, red = {}, {}
    for f in sorted(Path(folder).glob("pixel_check*/pixel_summary.csv")):
        for _, x in pd.read_csv(f).iterrows():
            pix[(x["target"], int(x["candidate"]))] = str(x["verdict"])
    for f in sorted(Path(folder).glob("reduction_check*/reduction_summary.csv")):
        for _, x in pd.read_csv(f).iterrows():
            red[(x["target"], int(x["candidate"]))] = str(x["verdict"])
    return {k: (pix.get(k, ""), red.get(k, "")) for k in set(pix) | set(red)}


def shortlist(folder):
    from enrich_candidates import period_matches
    folder = Path(folder)
    enriched, cands = folder / "enriched.csv", folder / "candidates_all.csv"
    if not cands.exists():
        cands = folder / "candidates.csv"
    if not enriched.exists():
        return None
    df = pd.read_csv(enriched, dtype={"flags": str, "enrich_flags": str})
    allrows = pd.read_csv(cands, dtype={"flags": str}) if cands.exists() else df
    notes = handled_notes()
    checks = later_checks(folder)
    out = []
    for _, r in df[df["status"] != "same_period"].iterrows():
        if flags_of(r.get("enrich_flags")) & set(SKIP_ENRICH):
            continue
        rel = allrows[(allrows["target"] == r["target"]) & (allrows["status"] == "same_period")
                      & allrows["flags"].fillna("").str.contains(f"related_to_{r['period_d']:.4f}",
                                                                regex=False)]
        why = signature(r, [x for _, x in rel.iterrows()])
        if not why:
            continue
        pix, red = checks.get((r["target"], int(r["candidate"])), ("", ""))
        if red == "reduction_dependent":                # not in other reductions: systematics
            continue
        cache = folder / "enrich_cache" / (re.sub(r"[^\w.-]+", "_", r["target"]) + ".json")
        info = json.loads(cache.read_text()) if cache.exists() else {}
        own = [v for v in info.get("variables", []) if v["sep"] < TARGET_ARCSEC]
        vsx = [v for v in own if v["cat"] == "VSX"]
        others = [v for v in own if v["cat"] != "VSX"]
        P = float(r["period_d"])

        def describe(v):
            p = f", P={v['period']:.5f} d" if np.isfinite(v["period"]) else ", no period"
            m = period_matches(P, [v])
            how = m[0][1].split(" = ", 1)[1].rsplit(",", 1)[0] if m else "no match"
            return f"{v['cat']} {v['name']} ({v['type'] or '?'}{p}{'' if p.endswith('period') else '; ' + how})"

        if not vsx:
            group = "new"
        else:
            confirmed = [v for v in vsx if np.isfinite(v["period"])
                         and any("alias" not in t for _, t in period_matches(P, [v]))]
            if confirmed:
                continue
            group = "period"
        # what we'd submit (an EB at this period) is already published there
        if any(v["cat"].startswith("T10k") and np.isfinite(v["period"])
               and any("alias" not in t for _, t in period_matches(P, [v])) for v in others):
            group = "published"
        blends = r.get("blend_suspects") if isinstance(r.get("blend_suspects"), str) else ""
        if "OFF-TARGET" in pix:                         # the dip may be on a neighbour
            plausible, faint = off_target_sources(pix, float(r["depth_ppm"]) * 1e-6, info.get("tmag"))
            if plausible:
                group = "blend"
                blends = "off target: " + "; ".join(plausible)
            else:                                       # no listed star could make the dip
                blends = ("pixel check off target, but no nearby star is bright enough to cause the dip"
                          + (f" ({faint})" if faint else ""))
        out.append(dict(group=group, star=r["target"], candidate=int(r["candidate"]), period_d=P,
                        depth_ppm=r["depth_ppm"], snr=r["snr"], signature=";".join(why),
                        vsx="; ".join(describe(v) for v in vsx),
                        other_catalogues="; ".join(describe(v) for v in others),
                        simbad_type=info.get("simbad", {}).get("otype", ""),
                        blend_suspects=blends,
                        note=notes.get(r["target"], ""), flags=r["flags"],
                        enrich_flags=r.get("enrich_flags") if isinstance(r.get("enrich_flags"), str) else ""))
    res = pd.DataFrame(out, columns=COLUMNS).sort_values(["group", "snr"], ascending=[True, False])
    res.to_csv(folder / "vsx_shortlist.csv", index=False)
    return res


def load_into_db(folder):
    """Replace the run's shortlist rows in tess.db, if the run has been ingested."""
    import tessdb
    folder = Path(folder).resolve()
    con = tessdb.connect()
    row = con.execute("SELECT run_id FROM runs WHERE folder = ?", (str(folder),)).fetchone()
    if not row:
        return False
    with con:
        con.execute("DELETE FROM vsx_shortlist WHERE run_id = ?", (row[0],))
        tessdb.load_vsx_shortlist(con, row[0], folder)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="results/<name> folder(s)")
    ap.add_argument("--no-db", action="store_true", help="only write the CSV")
    args = ap.parse_args()
    for f in args.runs:
        res = shortlist(f)
        if res is None:
            print(f"{f}: no enriched.csv (run enrich_candidates.py first)")
            continue
        loaded = not args.no_db and load_into_db(f)
        n = res["group"].value_counts().to_dict()
        print(f"{f}: {len(res)} for VSX review (new {n.get('new', 0)}, published {n.get('published', 0)}, period {n.get('period', 0)}, "
              f"blend {n.get('blend', 0)})"
              + (", loaded into tess.db" if loaded else ""))
        for _, r in res.iterrows():
            print(f"  {r['group']:6s} {r['star']:15s} P={r['period_d']:.5f} d  {r['depth_ppm']:>7.0f} ppm  "
                  f"SNR {r['snr']:5.1f}  {r['signature']}" + (f"  [{r['note'][:40]}]" if r["note"] else ""))


if __name__ == "__main__":
    main()
