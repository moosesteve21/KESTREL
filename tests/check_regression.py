#!/usr/bin/env python3
"""
tests/check_regression.py - compare a run of targets/regression.txt with the
expected outcomes in tests/regression_expect.csv.

Usage:
    .venv/bin/python tests/check_regression.py scratch/regression [--expect FILE]

Each expected signal is matched to a logged signal of the same star whose period is
within PERIOD_TOL (fractional). Expectations:
  clean       found, with status clean
  found       found, with any status except same_period
  not_clean   not clean: found with another status, or not logged at all
              (a rejected signal may be masked without a row)
plus, when given, a flag that must be among the signal's flags if it is found.
Prints one line per expectation and exits 1 if any fail.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PERIOD_TOL = 0.002


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="run folder (with candidates_all.csv)")
    ap.add_argument("--expect", default=str(HERE / "regression_expect.csv"))
    args = ap.parse_args()
    run = Path(args.run)
    cands = pd.read_csv(run / "candidates_all.csv")
    cands = cands[cands.status != "same_period"]
    done = {}
    if (run / "done_all.txt").exists():
        for line in (run / "done_all.txt").read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                done[parts[0]] = parts[1]
    exp = pd.read_csv(args.expect).fillna("")
    fails = 0
    for e in exp.itertuples():
        mine = cands[(cands.target == e.star) &
                     ((cands.period_d / e.period_d - 1).abs() < PERIOD_TOL)]
        got = mine.iloc[0] if len(mine) else None
        status = got["status"] if got is not None else None
        flags = str(got["flags"]) if got is not None and pd.notna(got["flags"]) else ""
        if e.star not in done:
            ok, why = False, "star not searched in this run"
        elif str(done[e.star]).startswith("error"):
            ok, why = False, f"search ended in {done[e.star]}"
        elif e.expect == "clean":
            ok, why = status == "clean", f"status {status}" if got is not None else "not found"
        elif e.expect == "found":
            ok, why = got is not None, f"status {status}" if got is not None else "not found"
        elif e.expect == "not_clean":
            ok, why = status != "clean", f"status {status}" if got is not None else "not logged"
        else:
            ok, why = False, f"unknown expectation {e.expect!r}"
        if ok and e.flag and got is not None and e.flag not in flags.split(";"):
            ok, why = False, f"{why}, missing flag {e.flag} (flags: {flags or 'none'})"
        fails += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {e.star:15s} P {e.period_d:<10g} expect {e.expect:9s} "
              f"{('+' + e.flag) if e.flag else '':28s} {why}   [{e.note}]")
    print(f"\n{len(exp) - fails}/{len(exp)} expectations met" + ("" if not fails else f"; {fails} FAILED"))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
