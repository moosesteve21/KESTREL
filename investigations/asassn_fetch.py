#!/usr/bin/env python3
"""Download precompiled ASAS-SN V-band light curves from the ASAS-SN Photometry Database
(https://asas-sn.osu.edu/photometry): no Sky Patrol computation needed.

Usage: .venv/bin/python investigations/asassn_fetch.py RA_DEG DEC_DEG OUTDIR [RADIUS_ARCMIN]

Searches the database around the position (default radius 0.5'), prints the sources found
(ID, position, separation, epochs, mean V, RMS, blend flag) and saves each one's light curve
as OUTDIR/asassn_<ID>.csv (columns hjd, camera, mag, mag_err, flux, flux_err). The database
holds V-band data, about 2012-2018 (DR9); the later g-band data need Sky Patrol v2.
Check that the listed position is your star before using the data.
"""
import html
import re
import sys
from pathlib import Path

import requests

BASE = "https://asas-sn.osu.edu/photometry"
HEADERS = {"User-Agent": "Mozilla/5.0 (PlanetFinder; asassn_fetch.py)"}
TIMEOUT = (15, 90)


def search(ra, dec, radius):
    r = requests.get(BASE, params={"ra": ra, "dec": dec, "radius": radius}, headers=HEADERS, timeout=TIMEOUT)
    r.raise_for_status()
    uuids = list(dict.fromkeys(re.findall(r'href="/photometry/([0-9a-f-]{36})"', r.text)))
    text = re.sub(r"<script.*?</script>|<style.*?</style>", " ", r.text, flags=re.S)
    cells = [c.strip() for c in html.unescape(re.sub(r"<[^>]+>", "\n", text)).splitlines() if c.strip()]
    rows = []
    for i, c in enumerate(cells):                 # ID, RA, Dec, distance, epochs, mean V, RMS, blend, source
        if re.fullmatch(r"AP\d+", c) and i + 8 < len(cells):
            rows.append(dict(zip(("id", "ra", "dec", "sep_arcsec", "epochs", "mean_v", "rms", "blend", "source"),
                                 cells[i:i + 9])))
    return list(zip(rows, uuids))


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    ra, dec, out = float(sys.argv[1]), float(sys.argv[2]), Path(sys.argv[3])
    radius = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    out.mkdir(parents=True, exist_ok=True)
    found = search(ra, dec, radius)
    if not found:
        print(f"no ASAS-SN photometry-database source within {radius}' of {ra} {dec}")
        return
    for row, uuid in found:
        print("  ".join(f"{k}={v}" for k, v in row.items()))
        r = requests.get(f"{BASE}/{uuid}.csv", headers=HEADERS, timeout=TIMEOUT)
        r.raise_for_status()
        path = out / f"asassn_{row['id']}.csv"
        path.write_text(r.text)
        print(f"  saved {path} ({len(r.text.splitlines()) - 1} rows)")


if __name__ == "__main__":
    main()
