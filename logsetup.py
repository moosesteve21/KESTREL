"""logsetup.py - one logging set-up for the follow-up scripts.

PLANETFINDER_LOG=debug (or info, warning) sets the level; the default is info, the
normal output. Debug adds the process details: each network query with its time and
retries, per-sector steps, the MCMC's convergence checks and TRICERATOPS's phases.

    PLANETFINDER_LOG=debug ./pipeline.sh targets/b0xx.txt results/b0xx 8   # a whole batch
    PLANETFINDER_LOG=debug .venv/bin/python pixel_check.py results/b0xx/enriched.csv
"""
import logging
import os
import time
from contextlib import contextmanager

FORMAT = "%(asctime)s %(levelname)s %(message)s"


def level():
    """The level PLANETFINDER_LOG asks for (default INFO)."""
    return getattr(logging, os.environ.get("PLANETFINDER_LOG", "info").upper(), logging.INFO)


def setup():
    logging.basicConfig(level=level(), format=FORMAT)
    for noisy in ("urllib3", "astroquery", "matplotlib", "PIL", "pyvo", "lightkurve"):  # their own debug output
        logging.getLogger(noisy).setLevel(max(level(), logging.INFO))


def debug_on():
    return level() <= logging.DEBUG


@contextmanager
def timed(what):
    """Log (debug) how long a block took: `with timed("TIC 123: TESScut search"): ...`"""
    t0 = time.time()
    try:
        yield
    finally:
        logging.debug(f"  [{time.time() - t0:6.1f} s] {what}")
