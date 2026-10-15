"""
vetting/transit_mcmc.py - physical transit fit (batman + emcee) for one candidate.
Runs in .venv-vet (Python 3.11); called by transit_fit.py, which prepares the data.

Usage: python vetting/transit_mcmc.py EVENTS.npz PARAMS.json
  EVENTS.npz: t, y, e, cad (days) - transits each normalised by a local baseline
  PARAMS.json: P, T0, depth, D (days), u1, u2 (limb-darkening prior means), steps
Prints one line of JSON with posterior medians and 16th/84th percentiles.

Model: circular orbit, quadratic limb darkening with Gaussian priors (width 0.1),
free epoch, period, Rp/R*, impact parameter b, log10 stellar density (uniform,
0.001-300 g/cm3) and a noise scale. Exposure smearing is modelled per cadence.
"""
import json
import sys
import warnings

warnings.filterwarnings("ignore")
import batman
import emcee
import numpy as np

G = 6.674e-8


CHUNK = 1000          # steps between convergence checks


def main():
    d = np.load(sys.argv[1])
    par = json.load(open(sys.argv[2]))
    t, y, e, cad = d["t"], d["y"], d["e"], d["cad"]
    P0, T00, D, depth0 = par["P"], par["T0"], par["D"], max(par["depth"], 1e-5)
    ul1, ul2 = par.get("u1", 0.3), par.get("u2", 0.3)
    max_steps = int(par.get("max_steps", par.get("steps", 30000)))
    groups = []
    for c in np.unique(np.round(cad * 1440, 1)):
        sel = np.abs(cad * 1440 - c) < 0.05
        pr = batman.TransitParams()
        pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.ecc, pr.w = T00, P0, np.sqrt(depth0), 15.0, 89.0, 0.0, 90.0
        pr.u, pr.limb_dark = [ul1, ul2], "quadratic"
        groups.append((sel, pr, batman.TransitModel(pr, t[sel], supersample_factor=5, exp_time=c / 1440)))
    k0 = np.sqrt(depth0)

    def model(th):
        t0, per, k, b, lrho, u1, u2, _ = th
        a = (G * 10 ** lrho * (per * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
        if b >= a:
            return None
        inc = np.degrees(np.arccos(b / a))
        out = np.empty_like(y)
        for sel, pr, m in groups:
            pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.u = t0, per, k, a, inc, [u1, u2]
            out[sel] = m.light_curve(pr)
        return out

    def log_prob(th):
        t0, per, k, b, lrho, u1, u2, lns = th
        if not (0.2 * k0 < k < 5 * k0 and k < 0.5 and 0 <= b < 1 + k and -3 < lrho < 2.5
                and 0 <= u1 <= 1 and 0 <= u2 <= 1 and u1 + u2 < 1 and abs(t0 - T00) < D
                and abs(per - P0) < 0.002 * P0 and -1 < lns < 1.5):
            return -np.inf
        m = model(th)
        if m is None:
            return -np.inf
        s = e * np.exp(lns)
        lp = -0.5 * ((u1 - ul1) / 0.1) ** 2 - 0.5 * ((u2 - ul2) / 0.1) ** 2
        return lp - 0.5 * np.sum(((y - m) / s) ** 2 + 2 * np.log(s))

    rng = np.random.default_rng(42)
    rho0 = par.get("rho_star") or 2.0
    p0 = np.array([T00, P0, k0, 0.3, np.log10(rho0), ul1, ul2, 0.0])
    nwalk = 32
    start = p0 + rng.normal(0, [D / 50, 1e-5 * P0, 0.05 * k0, 0.1, 0.1, 0.02, 0.02, 0.01], (nwalk, 8))
    # walkers spread over the whole impact-parameter range, grazing included, so the
    # sampler can reach a grazing solution (a V-shaped dip) instead of staying near b ~ 0.3
    start[:, 3] = rng.uniform(0.0, 1.0 + 0.8 * k0, nwalk)
    for i in range(nwalk):                       # move any walker the priors reject
        for _ in range(200):
            if np.isfinite(log_prob(start[i])):
                break
            start[i, 3] = rng.uniform(0.0, 1.0)
            start[i, 4] = np.log10(rho0) + rng.normal(0, 0.3)
    start[:, 5:7] = np.clip(start[:, 5:7], 0.01, 0.6)
    # Run until converged (emcee's recipe): the chain longer than 50 autocorrelation
    # times, and the estimate stable to 1%; checked every CHUNK steps, up to max_steps.
    # differential-evolution moves cope with the strong b / Rp/R* / density correlation
    # better than emcee's default stretch move (emcee docs, "Using different moves")
    sampler = emcee.EnsembleSampler(nwalk, 8, log_prob,
                                    moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    tau_old, converged, state = np.inf, False, start
    while sampler.iteration < max_steps:
        state = sampler.run_mcmc(state, CHUNK, progress=False)
        tau_all = sampler.get_autocorr_time(tol=0)
        tau = float(np.max(tau_all))
        if np.isfinite(tau) and sampler.iteration > 50 * tau and abs(tau_old - tau) / tau < 0.01:
            converged = True
            break
        tau_old = tau
    nstep = sampler.iteration
    burn = int(2 * tau) if np.isfinite(tau) else nstep // 2
    ch = sampler.get_chain(discard=min(burn, nstep // 2), thin=max(1, int(0.5 * tau)) if np.isfinite(tau) else 10,
                           flat=True)
    rho = 10 ** ch[:, 4]
    a = (G * rho * (ch[:, 1] * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    k, b = ch[:, 2], ch[:, 3]
    inc = np.arccos(np.clip(b / a, -1, 1))
    t14 = ch[:, 1] / np.pi * np.arcsin(np.clip(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) /
                                               (a * np.sin(inc)), 0, 1)) * 24

    def q(x):
        lo, mid, hi = np.percentile(x, [16, 50, 84])
        return [float(mid), float(lo), float(hi)]

    out = {"t0": q(ch[:, 0]), "period": q(ch[:, 1]), "rp_rs": q(k), "b": q(b), "log10_rho": q(ch[:, 4]),
           "rho": q(rho), "a_rs": q(a), "t14_h": q(t14), "depth_ppm": q(k ** 2 * 1e6),
           "u1": q(ch[:, 5]), "u2": q(ch[:, 6]), "noise_scale": q(np.exp(ch[:, 7])),
           "autocorr_steps": tau, "steps": nstep, "converged": converged, "samples": int(len(ch)),
           "n_points": int(len(t))}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
