"""Fouling scenario: how much of the energy survives a fouled membrane.

Extension to the report, not part of it. The report assumes clean membranes. In PRO the
effluent (feed) faces the porous support, so foulants from treated sewage can both
coat the membrane and clog the support. Two simple, time-averaged variants:

    "coating":           water permeability A x (1 - f)
    "coating + support": A x (1 - f) and structural parameter S x (1 + f)

for fouling severity f = 0 (clean) to 0.4. Membrane area is held at each reference design;
pressure and feed ratio are re-optimised (the operator adapts), as in s04_uncertainty.py.
The levels are illustrative, not measured: no pilot data for PRO on Gulf effluent was found
(bench forward-osmosis tests on treated wastewater show roughly 0-25 % flux decline).

A second part reruns the report's Monte Carlo (same 400 samples, same seed) with fouling
severity drawn from FOULING_RANGE and costs drawn from s06_economics.COSTS (about 5 minutes).

Writes results/fouling.csv, results/fouling_monte_carlo.csv, results/fouling_monte_carlo.json
and figures/fig12_fouling.png.
"""
import json
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd

from common import BASE_SITE, BASE_MACH, BLUE, ORANGE, INK2, RANGES, SEED, build, plt, save_json, savefig
from osmorec import CTA_HF, TFC_LAB, optimise, bulk_osmotic_difference, BAR
from s06_economics import COSTS, lcoe

FOULING_RANGE = (0.0, 0.15, 0.4)     # triangular (low, mode, high) for the combined Monte Carlo

MEMBRANES = {"cta": CTA_HF, "tfc": TFC_LAB}
LEVELS = [0.0, 0.1, 0.2, 0.3, 0.4]
VARIANTS = {"coating": lambda m, f: replace(m, A=m.A * (1 - f)),
            "coating + support": lambda m, f: replace(m, A=m.A * (1 - f), S=m.S * (1 + f))}


def worker(args):
    key, name, variant, f, sigma, x0 = args
    mem = VARIANTS[variant](MEMBRANES[key], f)
    d, b = optimise(BASE_SITE, mem, BASE_MACH, sigma, n_cells=60, x0=x0)
    return {"membrane": key, "design": name, "variant": variant, "fouling": f, "dP_bar": d.dP / BAR,
            "phi": d.phi, "net": b.net, "fraction": b.fraction_recovered, "area": b.area_per_product}


def mc_worker(args):
    """One Monte Carlo sample: the report's uncertain inputs plus fouling ("coating + support")."""
    key, params, f, sigma, x0 = args
    site, mem, mach = build(params, MEMBRANES[key])
    d, b = optimise(site, VARIANTS["coating + support"](mem, f), mach, sigma, n_cells=60, x0=x0)
    return {"fouling": f, "phi": d.phi, "net": b.net, "fraction": b.fraction_recovered, "area": b.area_per_product,
            "recovery": site.recovery}


if __name__ == "__main__":
    t0 = time.time()
    with open("results/designs.json") as fh:
        designs = json.load(fh)
    with open("results/design_study.json") as fh:
        bpp = json.load(fh)["site"]["brine_per_product"]
    dpi = bulk_osmotic_difference(BASE_SITE)
    jobs = []
    for key in MEMBRANES:
        for name in ("compact", "high_recovery"):
            r = designs[key][name]
            x0 = (r["dP_bar"] * BAR / dpi, r["phi"])
            jobs += [(key, name, v, f, r["sigma"], x0) for v in VARIANTS for f in LEVELS]
    with Pool() as pool:
        df = pd.DataFrame(pool.map(worker, jobs))

    clean = df[df.fouling == 0].set_index(["membrane", "design", "variant"]).fraction
    df["share_of_clean"] = [r.fraction / clean[(r.membrane, r.design, r.variant)] for r in df.itertuples()]
    base_cost = {p: b for p, (_, b, _, _) in COSTS.items()}
    df["lcoe_usd_kWh"] = [float(lcoe(r.net, r.area, r.phi * bpp, base_cost)["lcoe"]) for r in df.itertuples()]
    df.to_csv("results/fouling.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for ax, (key, name, title) in zip(axes, [("cta", "compact", "CTA hollow fibre, compact"),
                                             ("tfc", "high_recovery", "TFC lab membrane, high-recovery")]):
        for variant, color in zip(VARIANTS, (BLUE, ORANGE)):
            s = df[(df.membrane == key) & (df.design == name) & (df.variant == variant)]
            ax.plot(100 * s.fouling, 100 * s.fraction, marker="o", ms=4, color=color, label=variant)
        ax.set_title(title)
        ax.set_xlabel("Fouling severity f (%)")
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel("Plant electricity recovered (%)")
    axes[0].legend(title="Fouling acts on", fontsize=9, title_fontsize=9)
    fig.text(0.01, -0.04, "Illustrative fouling levels (time-averaged). Clean membrane at f = 0 matches the report.",
             color=INK2, fontsize=8)
    savefig(fig, "fig12_fouling.png")

    # Combined Monte Carlo: same physics samples as s04 (same seed), plus fouling and costs
    n_mc = 400
    summary, mc_parts = {}, []
    for key, name in (("cta", "compact"), ("tfc", "high_recovery")):
        r = designs[key][name]
        x0 = (r["dP_bar"] * BAR / dpi, r["phi"])
        rng = np.random.default_rng(SEED)
        samples = [{p: float(rng.triangular(lo, b0, hi)) for p, (lo, b0, hi, _) in RANGES.items()} for _ in range(n_mc)]
        f_s = np.random.default_rng(SEED + 7).triangular(*FOULING_RANGE, n_mc)
        with Pool() as pool:
            mc = pd.DataFrame(pool.map(mc_worker, [(key, s, float(f), r["sigma"], x0) for s, f in zip(samples, f_s)],
                                       chunksize=8))
        crng = np.random.default_rng(SEED + 8)
        c = {p: crng.triangular(lo, b0, hi, n_mc) for p, (lo, b0, hi, _) in COSTS.items()}
        bpp_s = (1 - mc.recovery) / mc.recovery
        mc["lcoe"] = pd.Series(lcoe(mc.net, mc.area, mc.phi * bpp_s, c)["lcoe"]).fillna(np.inf)
        mc["lcoe_no_pretreatment_capital"] = pd.Series(
            lcoe(mc.net, mc.area, mc.phi * bpp_s, c, with_pretreatment=False)["lcoe"]).fillna(np.inf)
        mc.insert(0, "design", name)
        mc.insert(0, "membrane", key)
        mc_parts.append(mc)
        q = lambda col: {k: float(v) for k, v in zip(("p05", "median", "p95"), mc[col].quantile([0.05, 0.5, 0.95]))}
        summary[f"{key}_{name}"] = {"fraction": q("fraction"), "lcoe": q("lcoe"),
                                    "lcoe_no_pretreatment_capital": q("lcoe_no_pretreatment_capital"),
                                    "prob_net_positive": float((mc.net > 0).mean()),
                                    "prob_lcoe_below_0.10": float((mc.lcoe < 0.10).mean())}
    pd.concat(mc_parts).to_csv("results/fouling_monte_carlo.csv", index=False)
    save_json("fouling_monte_carlo.json", {"n": n_mc, "fouling_variant": "coating + support",
                                           "fouling_triangular": FOULING_RANGE, "results": summary})

    pd.set_option("display.width", 200)
    print(json.dumps(summary, indent=1))
    print(df[["membrane", "design", "variant", "fouling", "fraction", "share_of_clean", "lcoe_usd_kWh"]]
          .round(3).to_string(index=False))
    print("done in %.0f s" % (time.time() - t0))
