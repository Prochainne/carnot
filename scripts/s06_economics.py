"""Economics: levelised cost of the electricity the PRO unit makes (US$ per kWh).

Extension to the report, not part of it: the report's numbers are energy only.
Reads results/designs.json and the Monte Carlo files from s04_uncertainty.py; runs no
new physics, so it takes seconds. Writes results/economics_designs.csv,
results/economics_monte_carlo.csv and results/economics.json.

Cost basis, per m3/day of desalinated product capacity:
    membrane capital  = area * membrane price / membrane share of PRO capital
                        (pressure vessels, turbine, pumps, pressure exchangers, piping)
    effluent capital  = effluent flow * pretreatment capital cost (UF or similar);
                        reported with and without it, since an existing reuse plant may
                        already treat the effluent well enough
    annual cost       = capital * (capital recovery factor + O&M share)
                        + membrane replacement (area * price / membrane life)
    annual energy     = net energy * 365 * availability
    LCOE              = annual cost / annual energy
Ranges and sources are in COSTS below and in data/SOURCES.md.
"""
import json

import numpy as np
import pandas as pd

from common import SEED, save_json, plt, savefig, BLUE, AQUA, YELLOW, ORANGE, INK2

# name: (low, base, high, label)
COSTS = {
    "membrane_usd_m2":   (5.0, 10.0, 30.0, "Membrane price (US$/m2)"),
    "membrane_share":    (0.5, 0.65, 0.8, "Membrane share of PRO capital"),
    "pretreat_usd_m3d":  (150.0, 400.0, 1000.0, "Effluent pretreatment capital (US$ per m3/d)"),
    "membrane_life_y":   (3.0, 5.0, 7.0, "Membrane life (years)"),
    "discount_rate":     (0.04, 0.07, 0.10, "Discount rate"),
    "om_share":          (0.02, 0.03, 0.04, "O&M per year (share of capital)"),
    "availability":      (0.85, 0.90, 0.95, "Availability"),
}
PLANT_LIFE_Y = 20
SOLAR_USD_KWH = 0.0135          # Al Dhafra PV tariff, Abu Dhabi, 2020
CHUNG_RO_BRINE_TSE = 0.44       # Chung et al. 2017, minimum LCOE for RO brine + wastewater


def crf(rate, years=PLANT_LIFE_Y):
    return rate * (1 + rate) ** years / ((1 + rate) ** years - 1)


def lcoe(net, area, effluent, c, with_pretreatment=True):
    """US$/kWh; all arguments per m3/day of product (net in kWh/m3). NaN if net <= 0.
    Works on scalars or numpy arrays."""
    net, area, effluent = np.asarray(net, float), np.asarray(area, float), np.asarray(effluent, float)
    cap_mem = area * c["membrane_usd_m2"] / c["membrane_share"]
    cap_pre = effluent * c["pretreat_usd_m3d"] if with_pretreatment else 0.0 * effluent
    k = crf(c["discount_rate"]) + c["om_share"]
    parts = {"membrane_and_equipment_capital": cap_mem * k,
             "membrane_replacement": area * c["membrane_usd_m2"] / c["membrane_life_y"],
             "pretreatment_capital": cap_pre * k}
    energy = net * 365.0 * c["availability"]
    with np.errstate(divide="ignore", invalid="ignore"):
        out = {p: np.where(net > 0, v / energy, np.nan) for p, v in parts.items()}
    out["lcoe"] = sum(out.values())
    return out


def breakeven_membrane_price(net, area, c, target):
    """Membrane price (US$/m2) at which LCOE without pretreatment capital equals target."""
    k = (crf(c["discount_rate"]) + c["om_share"]) / c["membrane_share"] + 1.0 / c["membrane_life_y"]
    return target * net * 365.0 * c["availability"] / (area * k)


if __name__ == "__main__":
    base = {p: b for p, (_, b, _, _) in COSTS.items()}
    with open("results/designs.json") as f:
        designs = json.load(f)
    with open("results/design_study.json") as f:
        bpp = json.load(f)["site"]["brine_per_product"]

    rows = []
    for key, d in designs.items():
        for name in ("compact", "high_recovery", "large_area"):
            r = d[name]
            eff = r["phi"] * bpp                               # m3 effluent per m3 product
            w, wo = lcoe(r["net"], r["area_per_product"], eff, base), lcoe(r["net"], r["area_per_product"], eff, base, False)
            rows.append({"membrane": key, "design": name, "fraction_recovered": r["fraction_recovered"],
                         "net_kWh_m3": r["net"], "area_m2_per_m3d": r["area_per_product"],
                         "effluent_m3_per_m3": eff,
                         **{f"usd_kWh_{p}": float(v) for p, v in w.items()},
                         "lcoe_no_pretreatment_capital": float(wo["lcoe"]),
                         "breakeven_membrane_usd_m2_solar": float(breakeven_membrane_price(
                             r["net"], r["area_per_product"], base, SOLAR_USD_KWH))})
    dfd = pd.DataFrame(rows)
    dfd.to_csv("results/economics_designs.csv", index=False)

    # Monte Carlo: physics samples from s04 paired with independently sampled costs
    rng = np.random.default_rng(SEED + 6)
    mc_rows, summary = [], {}
    for key, name in (("cta", "compact"), ("tfc", "high_recovery")):
        mc = pd.read_csv(f"results/monte_carlo_{key}_{name}.csv")
        n = len(mc)
        c = {p: rng.triangular(lo, b0, hi, n) for p, (lo, b0, hi, _) in COSTS.items()}
        sigma = designs[key][name]["sigma"]
        bpp_s = (1 - mc.recovery) / mc.recovery
        area = sigma * bpp_s / 24.0                            # m2 per (m3/day product)
        eff = mc.phi * bpp_s
        w, wo = lcoe(mc.net, area, eff, c), lcoe(mc.net, area, eff, c, False)
        part = pd.DataFrame({"membrane": key, "design": name, "net": mc.net, "fraction": mc.fraction, **c,
                             "lcoe": w["lcoe"], "lcoe_no_pretreatment_capital": wo["lcoe"]})
        mc_rows.append(part)
        s = {}
        for col in ("lcoe", "lcoe_no_pretreatment_capital"):
            v = part[col].fillna(np.inf)                        # net <= 0: never pays back
            q = v.quantile([0.05, 0.5, 0.95])
            s[col] = {"p05": q[0.05], "median": q[0.5], "p95": q[0.95],
                      "prob_below_solar": float((v < SOLAR_USD_KWH).mean()),
                      "prob_below_0.10": float((v < 0.10).mean())}
        s["share_net_not_positive"] = float(part.lcoe.isna().mean())
        summary[f"{key}_{name}"] = s
    pd.concat(mc_rows).to_csv("results/economics_monte_carlo.csv", index=False)

    save_json("economics.json", {"plant_life_years": PLANT_LIFE_Y, "base_costs": base,
                                 "ranges": {p: [lo, hi] for p, (lo, _, hi, _) in COSTS.items()},
                                 "comparison_usd_kWh": {"al_dhafra_solar_2020": SOLAR_USD_KWH,
                                                        "chung2017_ro_brine_wastewater": CHUNG_RO_BRINE_TSE},
                                 "monte_carlo": summary})
    # Figure 10: cost of each kWh by design, split into its parts
    parts = [("membrane_and_equipment_capital", "Membrane and equipment capital", BLUE),
             ("membrane_replacement", "Membrane replacement", AQUA),
             ("pretreatment_capital", "Effluent pretreatment capital", YELLOW)]
    names = {"cta": "CTA", "tfc": "TFC", "compact": "compact", "high_recovery": "high-recovery",
             "large_area": "large-area"}
    labels = [f"{names[r.membrane]} {names[r.design]} ({100 * r.fraction_recovered:.0f} %)" for r in dfd.itertuples()]
    fig, ax = plt.subplots(figsize=(8.5, 3.9))
    left = np.zeros(len(dfd))
    for col, label, color in parts:
        v = dfd[f"usd_kWh_{col}"].values
        ax.barh(labels, v, left=left, color=color, label=label, height=0.6)
        left += v
    ax.axvline(CHUNG_RO_BRINE_TSE, color=ORANGE, lw=1.2, ls="--")
    ax.text(CHUNG_RO_BRINE_TSE, 1.5, " Chung et al. 2017: 0.44", color=ORANGE, fontsize=8, va="center")
    ax.axvline(SOLAR_USD_KWH, color=INK2, lw=1.2)
    ax.text(SOLAR_USD_KWH, len(dfd) - 0.35, " Abu Dhabi solar: 0.0135", color=INK2, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Cost of the electricity produced (US$ per kWh), base-case costs")
    ax.set_title("What each kWh of osmotic power would cost")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", visible=False)
    savefig(fig, "fig10_economics.png")

    pd.set_option("display.width", 200)
    print(dfd[["membrane", "design", "fraction_recovered", "usd_kWh_lcoe", "lcoe_no_pretreatment_capital",
               "breakeven_membrane_usd_m2_solar"]].round(3).to_string(index=False))
    print(json.dumps(summary, indent=1, default=float))
