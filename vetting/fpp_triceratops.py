"""
vetting/fpp_triceratops.py - TRICERATOPS false-positive probabilities (Giacalone
et al. 2021) for one candidate. Runs in .venv-vet (Python 3.11); called by
fpp_check.py.

Usage: python vetting/fpp_triceratops.py PARAMS.json
  PARAMS.json: tic, P, T0 (BTJD), depth (fraction), events (npz from transit_fit.py),
  spoc_sectors, tess_spoc_sectors, draws, workdir
Prints one line of JSON: FPP, NFPP, scenario probabilities, sectors used.

Apertures: SPOC 2-min pipeline masks and TESS-SPOC optimal apertures (QLP publishes
none, so QLP-only sectors are left out of the flux-contribution step). The light
curve is the transit_fit.py events, folded and binned to 10 minutes. No contrast
curve. TRICERATOPS writes a Gaia query cache into workdir.
"""
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore")
import lightkurve as lk
import numpy as np
from astropy.io import fits


def main():
    par = json.load(open(sys.argv[1]))
    tic, P, T0, depth = int(par["tic"]), par["P"], par["T0"], par["depth"]
    os.makedirs(par["workdir"], exist_ok=True)
    os.chdir(par["workdir"])
    import triceratops.triceratops as tr
    apertures, sectors = [], []
    for s in par.get("spoc_sectors", []):
        try:
            tpf = lk.search_targetpixelfile(f"TIC {tic}", author="SPOC", exptime=120, sector=s).download()
            rows, cols = np.nonzero(tpf.pipeline_mask)
            if len(rows):
                apertures.append(np.column_stack([cols + tpf.column, rows + tpf.row]))
                sectors.append(s)
        except Exception:
            pass
    for s in par.get("tess_spoc_sectors", []):
        if s in sectors:
            continue
        try:
            lc = lk.search_lightcurve(f"TIC {tic}", author="TESS-SPOC", sector=s).download()
            with fits.open(lc.meta["FILENAME"]) as h:
                ap = h["APERTURE"]
                rows, cols = np.nonzero(ap.data & 2)
                if len(rows):
                    apertures.append(np.column_stack([cols + ap.header["CRVAL1P"], rows + ap.header["CRVAL2P"]]))
                    sectors.append(s)
        except Exception:
            pass
    if not sectors:
        print(json.dumps({"error": "no sector with a published aperture (SPOC or TESS-SPOC)"}))
        return
    target = tr.target(ID=tic, sectors=np.array(sectors))
    target.calc_depths(tdepth=depth, all_ap_pixels=apertures)
    d = np.load(par["events"])
    ph = (d["t"] - T0 + 0.5 * P) % P - 0.5 * P
    w = np.abs(ph) < 0.15
    edges = np.arange(-0.15, 0.15 + 1e-9, 10 / 1440)
    idx = np.digitize(ph[w], edges)
    keep = [i for i in range(1, len(edges)) if np.any(idx == i)]
    tb = np.array([ph[w][idx == i].mean() for i in keep])
    fb = np.array([d["y"][w][idx == i].mean() for i in keep])
    dur = par.get("t14_h", 2.0) / 24
    err = float(np.std(fb[np.abs(tb) > 0.75 * dur])) if np.sum(np.abs(tb) > 0.75 * dur) > 5 else float(np.std(fb))
    target.calc_probs(time=tb, flux_0=fb, flux_err_0=err, P_orb=P, N=int(par.get("draws", 1_000_000)),
                      parallel=True, verbose=0)
    probs = target.probs.sort_values("prob", ascending=False)
    nearby = target.stars[(target.stars["tdepth"] > 0) & (target.stars["tdepth"] <= 1) &
                          (target.stars["ID"].astype(str) != str(tic))]
    out = {"fpp": float(target.FPP), "nfpp": float(target.NFPP), "sectors": sectors,
           "scenarios": {f"{r.ID}:{r.scenario}": round(float(r.prob), 5) for r in probs.head(6).itertuples()},
           "prob_tp": float(probs[(probs.ID.astype(str) == str(tic)) & (probs.scenario == "TP")].prob.sum()),
           "n_nearby_possible": int(len(nearby)), "binned_points": int(len(tb)), "flux_err": err}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
