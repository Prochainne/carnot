"""Design study: optimum operating point versus membrane area, the two reference
designs per membrane, the energy waterfall, break-even in the feed-side energy
cost, the value of scarce effluent (Abu Dhabi scenarios) and the seawater-feed case.

Writes results/designs.json, results/design_study.json and several CSV files.
"""
import time
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from common import BASE_SITE, BASE_MACH, save_json
from osmorec import (CTA_HF, TFC_LAB, Membrane, Design, evaluate, optimise, optimise_pressure,
                     optimise_feed_value, ceiling, solve_module, bulk_osmotic_difference, BAR, KWH)

MEMBRANES = {"cta": CTA_HF, "tfc": TFC_LAB}
SIGMAS = np.unique(np.round(np.logspace(np.log10(8), np.log10(1500), 18), 1))
HIGH_RECOVERY_SHARE = 0.80      # "high-recovery" design: 80 % of the large-area net energy
SITE, MACH = BASE_SITE, BASE_MACH
DPI = bulk_osmotic_difference(SITE)


def row(d: Design, b):
    r = {"dP_bar": d.dP / BAR, "phi": d.phi, "sigma": d.sigma}
    r.update(b.as_dict())
    return r


def area_sweep(key):
    """Best (pressure, feed ratio) for each membrane area; then the two reference designs."""
    mem = MEMBRANES[key]
    rows, x0 = [], (0.45, 0.9)
    for s in SIGMAS[::-1]:                              # warm-start from the large-area end
        d, b = optimise(SITE, mem, MACH, float(s), x0=x0)
        x0 = (d.dP / DPI, d.phi)
        rows.append(row(d, b))
    df = pd.DataFrame(rows).sort_values("sigma").reset_index(drop=True)

    # compact design: maximum net electrical power per m2 of membrane
    i = int(df.power_density_net.idxmax())
    lo, hi = np.log(df.sigma[max(i - 1, 0)]), np.log(df.sigma[min(i + 1, len(df) - 1)])
    x0c = (df.dP_bar[i] * BAR / DPI, df.phi[i])
    r = minimize_scalar(lambda ls: -optimise(SITE, mem, MACH, float(np.exp(ls)), x0=x0c)[1].power_density_net,
                        bounds=(lo, hi), method="bounded", options={"xatol": 0.03})
    dc, bc = optimise(SITE, mem, MACH, float(np.exp(r.x)), x0=x0c)

    # high-recovery design: smallest area that yields 80 % of the large-area net energy
    target = HIGH_RECOVERY_SHARE * df.net.iloc[-1]
    sig_h = float(np.interp(target, df.net.values, df.sigma.values))
    j = int(np.argmin(np.abs(df.sigma.values - sig_h)))
    dh, bh = optimise(SITE, mem, MACH, sig_h, x0=(df.dP_bar[j] * BAR / DPI, df.phi[j]))
    return key, df, row(dc, bc), row(dh, bh)


def breakeven_curve(args):
    """Net energy versus the energy spent on each m3 of feed (pretreatment + conveyance)."""
    key, design_name, sigma, x0 = args
    mem, rows = MEMBRANES[key], []
    for e in np.round(np.arange(0.0, 0.725, 0.025), 3):
        d, b = optimise(SITE, mem, replace(MACH, e_pre_feed=float(e)), sigma, x0=x0)
        x0 = (d.dP / DPI, max(d.phi, 0.1))
        rows.append({"membrane": key, "design": design_name, "e_pre_feed": e, "fraction": b.fraction_recovered,
                     "net": b.net, "phi": d.phi, "dP_bar": d.dP / BAR, "net_per_m3_feed": b.net_per_m3_feed})
    return rows


def ideal_constant_pressure_limit(phi):
    """Hydraulic energy (kWh per m3 product) with an ideal membrane of unlimited area."""
    ideal = Membrane(A=1e-9, B=0.0, S=0.0, k=np.inf)

    def neg(f):
        r = solve_module(1.0, SITE.c_brine, phi, SITE.c_feed, f * DPI, 1e9, ideal, T=SITE.T, factor=SITE.salt_factor)
        return -f * DPI * r.permeate
    r = minimize_scalar(neg, bounds=(0.2, 0.9), method="bounded", options={"xatol": 1e-4})
    return -r.fun / KWH * SITE.brine_per_product, r.x * DPI / BAR


def as_design(r):
    return Design(r["dP_bar"] * BAR, r["phi"], r["sigma"])


if __name__ == "__main__":
    t0 = time.time()
    out = {"site": {"brine_salinity_g_kg": SITE.brine_salinity, "brine_per_product": SITE.brine_per_product,
                    "osmotic_difference_bar": DPI / BAR, "c_brine_mol_m3": SITE.c_brine,
                    "c_feed_mol_m3": SITE.c_feed}}
    designs = {}
    with Pool() as pool:
        for key, df, compact, high in pool.map(area_sweep, list(MEMBRANES)):
            df.to_csv(f"results/area_sweep_{key}.csv", index=False)
            designs[key] = {"membrane": MEMBRANES[key].name, "compact": compact, "high_recovery": high,
                            "large_area": df.iloc[-1].to_dict(),
                            "max_gross_power_density": float(df.power_density_gross.max())}
        save_json("designs.json", designs)
        jobs = [(k, n, designs[k][n]["sigma"], (designs[k][n]["dP_bar"] * BAR / DPI, designs[k][n]["phi"]))
                for k in MEMBRANES for n in ("compact", "high_recovery")]
        be = pd.DataFrame([r for rows in pool.map(breakeven_curve, jobs) for r in rows])
    be.to_csv("results/breakeven_pretreatment.csv", index=False)
    out["breakeven"] = {}
    for (k, n), f in be.groupby(["membrane", "design"]):
        out["breakeven"][f"{k}_{n}"] = float(np.interp(0.0, f.net.values[::-1], f.e_pre_feed.values[::-1]))

    # energy waterfall, from the thermodynamic ceiling down to net electricity
    out["waterfall"] = {}
    for k in MEMBRANES:
        for n in ("compact", "high_recovery"):
            r = designs[k][n]
            ceil = ceiling(SITE, r["phi"])
            ideal_cp, dp_ideal = ideal_constant_pressure_limit(r["phi"])
            out["waterfall"][f"{k}_{n}"] = {
                "reversible_ceiling": ceil, "ideal_constant_pressure": ideal_cp,
                "ideal_constant_pressure_dP_bar": dp_ideal,
                "hydraulic_at_design": r["gross_turbine"] / MACH.eta_turbine,
                "after_turbine": r["gross_turbine"], "loss_px": r["loss_px"],
                "loss_draw_pump": r["loss_draw_pump"], "loss_feed_pump": r["loss_feed_pump"],
                "loss_pretreatment": r["loss_pretreatment"], "net": r["net"],
                "net_share_of_ceiling": r["net"] / ceil}

    # concentration and flux profiles along the module (CTA, compact design)
    d = as_design(designs["cta"]["compact"])
    res = solve_module(1.0, SITE.c_brine, d.phi, SITE.c_feed, d.dP, d.sigma * 3600, CTA_HF,
                       T=SITE.T, factor=SITE.salt_factor, n_cells=200)
    mid = 0.5 * (res.area_profile[:-1] + res.area_profile[1:]) / res.area_profile[-1]
    pd.DataFrame({"area_fraction": mid, "cD_mol_m3": res.cD, "cF_mol_m3": res.cF,
                  "Jw_LMH": res.Jw * 3.6e6}).to_csv("results/module_profile_cta.csv", index=False)

    # scarce effluent: best energy per m3 of effluent, and Abu Dhabi scenarios
    PRODUCT, EMIRATE_TSE = 909200.0, 1.0e6               # m3/day
    brine_day = PRODUCT * SITE.brine_per_product
    ro_mw = SITE.sec_ro * PRODUCT / 24.0 / 1e3
    out["abu_dhabi"] = {"product_m3_day": PRODUCT, "brine_m3_day": brine_day, "ro_power_MW": ro_mw,
                        "emirate_effluent_m3_day": EMIRATE_TSE}
    out["effluent_value"], rows = {}, []
    for k, mem in MEMBRANES.items():
        for n in ("compact", "high_recovery"):
            r = designs[k][n]
            d0, b0 = as_design(r), evaluate(SITE, mem, MACH, as_design(r))
            sigma_feed = r["sigma"] / r["phi"]           # membrane area per m3/h of effluent
            dv, bv = optimise_feed_value(SITE, mem, MACH, sigma_feed, x0=(r["dP_bar"] * BAR / DPI, r["phi"] * 0.6))
            out["effluent_value"][f"{k}_{n}"] = {
                "sigma_feed": sigma_feed, "best_kWh_per_m3_effluent": bv.net_per_m3_feed, "phi": dv.phi,
                "dP_bar": dv.dP / BAR, "at_design_kWh_per_m3_effluent": b0.net_per_m3_feed,
                "share_of_desalinating_1_m3": bv.net_per_m3_feed / SITE.sec_ro}
            for label, tse in [("Surplus effluent only (20 % of the emirate's 1.0 Mm3/d)", 0.2 * EMIRATE_TSE),
                               ("Half of the emirate's effluent", 0.5 * EMIRATE_TSE),
                               ("All of the emirate's effluent", EMIRATE_TSE),
                               ("Unconstrained (design optimum)", r["phi"] * brine_day)]:
                phi_avail = tse / brine_day
                if phi_avail >= r["phi"] * 0.999:        # effluent is not limiting
                    bb, used, dp = b0, r["phi"] * brine_day, d0.dP
                elif phi_avail <= dv.phi:                # use part of the brine at the best ratio
                    bb, used, dp = bv, tse, dv.dP
                else:                                    # all the brine, all the available effluent
                    dd, bb = optimise_pressure(SITE, mem, MACH, phi_avail, sigma_feed * phi_avail)
                    used, dp = tse, dd.dP
                net_mw = bb.net_per_m3_feed * used / 24.0 / 1e3
                rows.append({"membrane": k, "design": n, "scenario": label, "effluent_m3_day": used,
                             "share_of_emirate_effluent": used / EMIRATE_TSE,
                             "kWh_per_m3_effluent": bb.net_per_m3_feed, "net_MW": net_mw,
                             "fraction_of_ro_demand": net_mw / ro_mw,
                             "membrane_area_million_m2": sigma_feed * used / 24.0 / 1e6, "dP_bar": dp / BAR})
    pd.DataFrame(rows).to_csv("results/abu_dhabi_scenarios.csv", index=False)

    # seawater instead of effluent as the feed
    sw_site = replace(SITE, feed_salinity=SITE.seawater_salinity)
    rows = []
    for k, mem in MEMBRANES.items():
        for n in ("compact", "high_recovery"):
            for e_pre in (0.0, 0.1):
                for phi in (0.25, 0.5, 1.0, 2.0):
                    dd, bb = optimise_pressure(sw_site, mem, replace(MACH, e_pre_feed=e_pre), phi, designs[k][n]["sigma"])
                    rows.append({"membrane": k, "design": n, "e_pre_feed": e_pre, "phi": phi, "dP_bar": dd.dP / BAR,
                                 "gross": bb.gross_turbine, "net": bb.net, "fraction": bb.fraction_recovered})
    sw = pd.DataFrame(rows)
    sw.to_csv("results/seawater_feed.csv", index=False)
    out["seawater_feed"] = {"osmotic_difference_bar": bulk_osmotic_difference(sw_site) / BAR,
                            "best_fraction": float(sw.fraction.max()),
                            "best_fraction_no_pretreatment": float(sw[sw.e_pre_feed == 0].fraction.max())}
    save_json("design_study.json", out)
    for k in designs:
        for n in ("compact", "high_recovery", "large_area"):
            r = designs[k][n]
            print(k, n, {x: round(r[x], 3) for x in ("sigma", "dP_bar", "phi", "net", "fraction_recovered",
                  "power_density_gross", "power_density_net", "area_per_product", "discharge_salinity", "net_per_m3_feed")})
    print("done in %.0f s" % (time.time() - t0))
