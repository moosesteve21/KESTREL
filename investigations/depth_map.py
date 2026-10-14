#!/usr/bin/env python3
"""Where on the sky is a candidate's dip? A per-source pixel depth map from TESScut
cutouts; the method is described in depth_map.py (repo root), which the pipeline's
depth_map_check.py also uses.

Usage: .venv/bin/python investigations/depth_map.py TIC P_DAYS T0_BTJD T14_HOURS OUTDIR
TIC is the number alone ("TIC " is accepted too). Writes depth_map.png and depth_map.txt
in OUTDIR. DEPTHMAP_SKIP=63,64 leaves sectors out (for example one dominated by
scattered light).
"""
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import depth_map as dm  # noqa: E402

if len(sys.argv) != 6:
    sys.exit(__doc__)
TIC = int(re.sub(r"\D", "", sys.argv[1]))
P, T0, T14 = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) / 24
OUT = Path(sys.argv[5]); OUT.mkdir(parents=True, exist_ok=True)
SKIP = {int(x) for x in os.environ.get("DEPTHMAP_SKIP", "").split(",") if x.strip()}
lines = []


def say(s=""):
    print(s)
    lines.append(s)


res = dm.depth_map(TIC, P, T0, T14, skip=SKIP, say=say)
say()
dm.ranking(res, say=say)
(OUT / "depth_map.txt").write_text("\n".join(lines) + "\n")
dm.plot(res, OUT / "depth_map.png")
