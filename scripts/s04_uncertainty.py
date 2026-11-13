"""Uncertainty: one-at-a-time sensitivity (tornado), Monte Carlo, and the break-even map.

Usage:
    python scripts/s04_uncertainty.py <cta|tfc> <compact|high_recovery> [n_samples]
    python scripts/s04_uncertainty.py map

The membrane area is held at the reference design; the applied pressure and the
feed ratio are re-optimised for every parameter set (the operator adapts).
Requires results/designs.json from s03_design.py.
"""
import json
import sys
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import BASE_SITE, BASE_MACH, RANGES, SEED, build, save_json
from osmorec import CTA_HF, TFC_LAB, optimise, bulk_osmotic_difference, BAR

MEMBRANES = {"cta": CTA_HF, "tfc": TFC_LAB}


def worker(args):
    key, params, sigma, x0 = args
    site, mem, mach = build(params, MEMBRANES[key])
    d, b = optimise(site, mem, mach, sigma, n_cells=60, x0=x0)
    return {**params, "dP_bar": d.dP / BAR, "phi": d.phi, "net": b.net, "fraction": b.fraction_recovered,
            "net_per_m3_feed": b.net_per_m3_feed, "pd_gross": b.power_density_gross}


def map_worker(args):
    key, e_pre, s_feed, sigma, x0 = args
    site = replace(BASE_SITE, feed_salinity=s_feed)
    d, b = optimise(site, MEMBRANES[key], replace(BASE_MACH, e_pre_feed=e_pre), sigma, n_cells=60, x0=x0)
    return {"e_pre_feed": e_pre, "feed_salinity": s_feed, "fraction": b.fraction_recovered, "net": b.net,
            "phi": d.phi, "dP_bar": d.dP / BAR}


if __name__ == "__main__":
    t0 = time.time()
    with open("results/designs.json") as f:
        designs = json.load(f)
    dpi = bulk_osmotic_difference(BASE_SITE)
    if sys.argv[1] == "map":
        r = designs["cta"]["compact"]
        x0 = (r["dP_bar"] * BAR / dpi, r["phi"])
        e_grid = np.linspace(0.0, 0.6, 13)
        s_grid = [0.5, 1, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 42]
        jobs = [("cta", float(e), float(s), r["sigma"], x0) for s in s_grid for e in e_grid]
        with Pool() as pool:
            pd.DataFrame(pool.map(map_worker, jobs, chunksize=4)).to_csv("results/breakeven_map_cta_compact.csv", index=False)
    else:
        key, name = sys.argv[1], sys.argv[2]
        n_mc = int(sys.argv[3]) if len(sys.argv) > 3 else 400
        r = designs[key][name]
        sigma, x0, base = r["sigma"], (r["dP_bar"] * BAR / dpi, r["phi"]), r["fraction_recovered"]
        rng = np.random.default_rng(SEED)
        samples = [{p: float(rng.triangular(lo, b0, hi)) for p, (lo, b0, hi, _) in RANGES.items()} for _ in range(n_mc)]
        jobs = []
        for p, (lo, _, hi, _) in RANGES.items():
            jobs += [(key, {p: lo}, sigma, x0), (key, {p: hi}, sigma, x0)]
        with Pool() as pool:
            res_t = pool.map(worker, jobs, chunksize=2)
            mc = pd.DataFrame(pool.map(worker, [(key, s, sigma, x0) for s in samples], chunksize=8))
        rows = [{"parameter": p, "label": label, "low": lo, "base": b0, "high": hi,
                 "fraction_low": res_t[2 * i]["fraction"], "fraction_high": res_t[2 * i + 1]["fraction"],
                 "fraction_base": base}
                for i, (p, (lo, b0, hi, label)) in enumerate(RANGES.items())]
        td = pd.DataFrame(rows)
        td["swing"] = (td.fraction_high - td.fraction_low).abs()
        td.sort_values("swing", ascending=False).to_csv(f"results/tornado_{key}_{name}.csv", index=False)
        mc.to_csv(f"results/monte_carlo_{key}_{name}.csv", index=False)
        q = mc.fraction.quantile([0.05, 0.25, 0.5, 0.75, 0.95])
        rho = {p: float(spearmanr(mc[p], mc.fraction)[0]) for p in RANGES}
        save_json(f"uncertainty_{key}_{name}.json", {
            "membrane": key, "design": name, "sigma": sigma, "n": n_mc, "seed": SEED, "fraction_base": base,
            "fraction_p05": q[0.05], "fraction_p25": q[0.25], "fraction_median": q[0.5],
            "fraction_p75": q[0.75], "fraction_p95": q[0.95], "fraction_mean": mc.fraction.mean(),
            "fraction_min": mc.fraction.min(), "fraction_max": mc.fraction.max(),
            "prob_net_positive": float((mc.net > 0).mean()),
            "prob_fraction_above_5pct": float((mc.fraction > 0.05).mean()),
            "prob_fraction_above_10pct": float((mc.fraction > 0.10).mean()),
            "spearman": dict(sorted(rho.items(), key=lambda kv: -abs(kv[1])))})
    print("done in %.0f s" % (time.time() - t0))
