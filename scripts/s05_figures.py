"""Draw the result figures from the files in results/ (no model runs here)."""
import json

import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from common import (savefig, plt, BLUE, ORANGE, AQUA, INK, INK2, MUTED, GRID, AXIS, SURFACE, RANGES)

with open("results/designs.json") as f:
    DES = json.load(f)
with open("results/design_study.json") as f:
    STUDY = json.load(f)
NAMES = {"cta": "Commercial CTA hollow fibre", "tfc": "Lab thin-film composite"}
COL = {"cta": BLUE, "tfc": AQUA}


# --- system schematic -----------------------------------------------------------
def schematic():
    fig, ax = plt.subplots(figsize=(10.5, 4.7))
    ax.set_xlim(0, 122)
    ax.set_ylim(2, 50)
    ax.axis("off")

    def box(x, y, w, h, text, fc="#ffffff", ec=AXIS, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=fc, ec=ec, lw=1.4))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9.5,
                fontweight="bold" if bold else "normal", color=INK)

    def arrow(p, q, color=INK2, lw=1.8):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=13, lw=lw, color=color,
                                     shrinkA=0, shrinkB=0))

    def note(x, y, text, ha="center", va="bottom", color=INK2):
        ax.text(x, y, text, ha=ha, va=va, fontsize=8.5, color=color)

    box(14, 30, 22, 10, "Seawater RO plant\n3.0 kWh per m³", bold=True)
    arrow((1, 35), (14, 35))
    note(7, 36, "Seawater\n42 g/kg")
    arrow((25, 40), (25, 48))
    note(27, 44, "Fresh water, 1 m³", ha="left", va="center")
    arrow((36, 35), (50, 35), color=ORANGE)
    note(43, 36, "Brine 1.38 m³\n72 g/kg")
    box(50, 30, 16, 10, "Pressure\nexchanger")
    arrow((66, 35), (78, 35), color=ORANGE, lw=2.6)
    note(72, 36, "brine at\npressure ΔP")
    ax.add_patch(FancyBboxPatch((78, 12), 28, 28, boxstyle="round,pad=0.3,rounding_size=1.2", fc="#ffffff",
                                ec=AXIS, lw=1.4))
    ax.plot([78, 106], [26, 26], color=INK, lw=1.2, ls=(0, (4, 2)))
    ax.text(92, 36.5, "PRO membrane module", ha="center", va="center", fontsize=9.5, fontweight="bold")
    note(92, 31, "brine side, high pressure", va="center")
    note(92, 15, "effluent side, low pressure", va="center")
    for x in (82, 102):
        arrow((x, 19.5), (x, 29.5), color=BLUE, lw=1.6)
    note(92, 22.3, "water crosses the\nmembrane by osmosis", va="center", color=BLUE)
    arrow((120, 19), (106, 19), color=BLUE)
    note(113.5, 20, "Treated sewage\neffluent, ~1 g/kg")
    arrow((92, 12), (92, 4), color=BLUE)
    note(93.5, 7.5, "concentrated effluent out", ha="left", va="center")
    ax.plot([106, 113, 113, 58], [35, 35, 46.5, 46.5], color=ORANGE, lw=2.6, solid_joinstyle="round")
    arrow((58, 46.5), (58, 40), color=ORANGE, lw=2.6)
    note(88, 47.3, "diluted brine, still at pressure ΔP")
    arrow((58, 30), (58, 22), color=ORANGE)
    note(56.5, 26, "extra volume (= the\nwater that crossed)", ha="right", va="center")
    box(48, 12, 20, 10, "Turbine and\ngenerator", fc="#fff3ec", ec=ORANGE, bold=True)
    ax.plot([48, 25], [17, 17], color=INK, lw=1.8)
    arrow((25, 17), (25, 30), color=INK)
    note(36.5, 18, "electricity")
    arrow((58, 12), (58, 4), color=ORANGE)
    note(59.5, 7.5, "diluted brine to sea", ha="left", va="center")
    savefig(fig, "fig0_schematic.png")


# --- trade-off between recovered energy and membrane area -----------------------
def tradeoff():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for key in ("cta", "tfc"):
        df = pd.read_csv(f"results/area_sweep_{key}.csv")
        df = df[df.net > -0.02]
        axes[0].plot(df.area_per_product, 100 * df.fraction_recovered, color=COL[key], label=NAMES[key])
        axes[1].plot(df.area_per_product, df.power_density_gross, color=COL[key], label=NAMES[key])
        for name, mk in (("compact", "o"), ("high_recovery", "D")):
            r = DES[key][name]
            for ax, y in ((axes[0], 100 * r["fraction_recovered"]), (axes[1], r["power_density_gross"])):
                ax.plot(r["area_per_product"], y, mk, color=COL[key], ms=8, mec="white", mew=1.5, zorder=5)
        axes[0].axhline(100 * DES[key]["large_area"]["fraction_recovered"], color=COL[key], lw=1, ls=":")
    r = DES["cta"]["compact"]
    axes[0].annotate("compact design\n(most power per m²)", (r["area_per_product"], 100 * r["fraction_recovered"]),
                     xytext=(9, 1.0), fontsize=8.5, color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED))
    r = DES["cta"]["high_recovery"]
    axes[0].annotate("high-recovery design\n(80 % of the limit)", (r["area_per_product"], 100 * r["fraction_recovered"]),
                     xytext=(13, 17.0), fontsize=8.5, ha="right", color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED))
    axes[0].text(0.56, 100 * DES["cta"]["large_area"]["fraction_recovered"] + 0.3, "unlimited-area limits (dotted)",
                 ha="left", fontsize=8.5, color=INK2)
    axes[1].axhline(5, color=ORANGE, lw=1.5)
    axes[1].text(85, 5.3, "5 W/m² viability threshold cited in the literature", ha="right", fontsize=8.5, color=INK2)
    axes[0].set_ylabel("RO electricity recovered (%)")
    axes[1].set_ylabel("Gross membrane power density (W/m²)")
    axes[0].set_title("More membrane recovers more energy ...")
    axes[1].set_title("... but each square metre does less work")
    axes[0].axhline(0, color=AXIS, lw=1)
    axes[0].set_ylim(-1, 21)
    axes[1].set_ylim(0, 12)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlim(0.5, 90)
        ax.set_xlabel("Membrane area per unit of RO capacity (m² per m³/day)")
    axes[1].legend(loc="upper right", fontsize=8.5)
    axes[0].text(0.56, 15.2, "circle: compact design\ndiamond: high-recovery design", fontsize=8.5, color=INK2)
    fig.tight_layout(w_pad=2.5)
    savefig(fig, "fig4_area_tradeoff.png")


# --- waterfall: from the thermodynamic ceiling to net electricity ---------------
def waterfall():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), sharey=True, sharex=True)
    for ax, (key, name, title) in zip(axes, (("cta", "compact", "Commercial membrane, compact design"),
                                             ("cta", "high_recovery", "Commercial membrane, high-recovery design"))):
        w = STUDY["waterfall"][f"{key}_{name}"]
        steps = [("Reversible ceiling", w["reversible_ceiling"], None),
                 ("Constant pressure", None, w["reversible_ceiling"] - w["ideal_constant_pressure"]),
                 ("Real membrane, finite area", None, w["ideal_constant_pressure"] - w["hydraulic_at_design"]),
                 ("Turbine", None, w["hydraulic_at_design"] - w["after_turbine"]),
                 ("Pressure exchanger", None, w["loss_px"]),
                 ("Pumping", None, w["loss_draw_pump"] + w["loss_feed_pump"]),
                 ("Effluent pretreatment", None, w["loss_pretreatment"]),
                 ("Net electricity", w["net"], None)]
        level = 0.0
        for i, (lab, total, loss) in enumerate(steps):
            y = len(steps) - 1 - i
            if total is not None:
                ax.barh(y, total, color=BLUE, height=0.62)
                ax.text(total + 0.02, y, f"{total:.2f}", va="center", fontsize=8.5, color=INK)
                level = total
            else:
                ax.barh(y, loss, left=level - loss, color=ORANGE, height=0.62)
                ax.text(level + 0.02, y, f"−{loss:.2f}", va="center", fontsize=8.5, color=INK)
                level -= loss
        ax.set_yticks(range(len(steps)))
        ax.set_yticklabels([s[0] for s in steps][::-1], fontsize=9)
        ax.set_title(title)
        ax.grid(axis="y", visible=False)
        ax.set_axisbelow(True)
        ax.set_xlabel("Energy (kWh per m³ of desalinated water)")
        ax.set_xlim(0, 1.6)
    fig.tight_layout(w_pad=2.0)
    savefig(fig, "fig5_waterfall.png")


# --- profiles along the module ---------------------------------------------------
def profiles():
    p = pd.read_csv("results/module_profile_cta.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6))
    axes[0].plot(p.area_fraction, p.cD_mol_m3 * 0.058443, color=ORANGE, label="Brine (draw)")
    axes[0].plot(p.area_fraction, p.cF_mol_m3 * 0.058443, color=BLUE, label="Effluent (feed)")
    axes[0].set_ylabel("Bulk salt concentration (g/L)")
    axes[0].set_title("Brine dilutes, effluent concentrates")
    axes[0].legend(loc="center right", fontsize=8.5)
    axes[0].annotate("", (0.98, 3), (0.72, 3), arrowprops=dict(arrowstyle="<|-", color=BLUE))
    axes[0].text(0.85, 5, "effluent flow", ha="center", fontsize=8.5, color=INK2)
    axes[0].annotate("", (0.30, 46), (0.04, 46), arrowprops=dict(arrowstyle="-|>", color=ORANGE))
    axes[0].text(0.17, 48, "brine flow", ha="center", fontsize=8.5, color=INK2)
    axes[1].plot(p.area_fraction, p.Jw_LMH, color=AQUA)
    axes[1].set_ylabel("Water flux (L m$^{-2}$ h$^{-1}$)")
    axes[1].set_title("Local water flux")
    for ax in axes:
        ax.set_xlabel("Position along the membrane (0 = brine inlet)")
        ax.set_ylim(bottom=0)
        ax.set_xlim(0, 1)
    savefig(fig, "fig6_module_profiles.png")


# --- break-even --------------------------------------------------------------------
def breakeven():
    be = pd.read_csv("results/breakeven_pretreatment.csv")
    mp = pd.read_csv("results/breakeven_map_cta_compact.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    ax = axes[0]
    for key in ("cta", "tfc"):
        for name, ls in (("compact", "-"), ("high_recovery", "--")):
            f = be[(be.membrane == key) & (be.design == name) & (be.phi > 0.05)]
            ax.plot(f.e_pre_feed, 100 * f.fraction, color=COL[key], ls=ls,
                    label=f"{NAMES[key]}, {'compact' if name == 'compact' else 'high recovery'}")
    ax.axhline(0, color=INK, lw=1)
    ax.axvspan(0.05, 0.30, color=MUTED, alpha=0.12, lw=0)
    ax.text(0.175, -2.6, "range assumed", ha="center", fontsize=8.5, color=INK2)
    ax.set_xlabel("Energy spent on each m³ of effluent (kWh/m³)")
    ax.set_ylabel("RO electricity recovered (%)")
    ax.set_title("Net gain versus the energy cost of the feed")
    ax.set_ylim(-3, 23)
    ax.legend(loc="upper right", fontsize=7.8)
    ax = axes[1]
    piv = mp.pivot(index="feed_salinity", columns="e_pre_feed", values="fraction") * 100
    X, Y = np.meshgrid(piv.columns.values, piv.index.values)
    cs = ax.contourf(X, Y, np.clip(piv.values, 0, None), levels=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                     cmap=plt.cm.Blues, extend="neither")
    ax.contourf(X, Y, piv.values, levels=[-100, 0], colors=["#f0efec"])
    ax.contour(X, Y, piv.values, levels=[0], colors=[INK], linewidths=2)
    cb = fig.colorbar(cs, ax=ax, pad=0.02)
    cb.set_label("RO electricity recovered (%)")
    ax.set_yscale("log")
    ax.set_yticks([0.5, 1, 2, 5, 10, 20, 42])
    ax.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "42"])
    ax.text(0.52, 24, "no net gain", fontsize=9, color=INK2, ha="center")
    ax.text(0.59, 36, "seawater = 42", fontsize=8.5, color=INK2, ha="right", va="center")
    ax.set_xlabel("Energy spent on each m³ of feed (kWh/m³)")
    ax.set_ylabel("Feed salinity (g/kg)")
    ax.set_title("Break-even map (commercial, compact)")
    ax.grid(False)
    fig.tight_layout(w_pad=2.0)
    savefig(fig, "fig7_breakeven.png")


# --- tornado ------------------------------------------------------------------------
def tornado():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2))
    for ax, (key, name, title) in zip(axes, (("cta", "compact", "Commercial membrane, compact design"),
                                             ("tfc", "high_recovery", "Lab membrane, high-recovery design"))):
        t = pd.read_csv(f"results/tornado_{key}_{name}.csv").head(12).iloc[::-1]
        base = 100 * t.fraction_base.iloc[0]
        y = np.arange(len(t))
        lo = 100 * t[["fraction_low", "fraction_high"]].min(axis=1).values
        hi = 100 * t[["fraction_low", "fraction_high"]].max(axis=1).values
        ax.barh(y, hi - base, left=base, color=BLUE, height=0.6)
        ax.barh(y, lo - base, left=base, color=ORANGE, height=0.6)
        ax.axvline(base, color=INK, lw=1)
        ax.set_yticks(y)
        ax.set_yticklabels(t.label.values, fontsize=8.3)
        ax.set_xlabel("RO electricity recovered (%)")
        ax.set_title(f"{title} (base {base:.1f} %)")
        ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.text(0.5, -0.02, "Each input moved alone across its range. Blue: better than the base case. Orange: worse.",
             ha="center", fontsize=8.8, color=INK2)
    savefig(fig, "fig8_tornado.png")


# --- Monte Carlo -------------------------------------------------------------------
def monte_carlo():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    for ax, (key, name, title) in zip(axes, (("cta", "compact", "Commercial membrane, compact design"),
                                             ("tfc", "high_recovery", "Lab membrane, high-recovery design"))):
        mc = pd.read_csv(f"results/monte_carlo_{key}_{name}.csv")
        with open(f"results/uncertainty_{key}_{name}.json") as f:
            u = json.load(f)
        ax.hist(100 * mc.fraction, bins=28, color=COL[key], edgecolor=SURFACE, linewidth=1)
        for q, lab in ((u["fraction_p05"], "5 %"), (u["fraction_median"], "median"), (u["fraction_p95"], "95 %")):
            ax.axvline(100 * q, color=INK, lw=1, ls="--")
            ax.text(100 * q, ax.get_ylim()[1] * 1.0, f"{lab}\n{100 * q:.1f}", ha="center", va="bottom", fontsize=8.3, color=INK2)
        ax.axvline(100 * u["fraction_base"], color=ORANGE, lw=2)
        ax.set_xlabel("RO electricity recovered (%)")
        ax.set_title(title, pad=26)
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel(f"Samples (of {len(mc)})")
    axes[0].text(0.02, 0.9, "orange line: base case", transform=axes[0].transAxes, fontsize=8.3, color=INK2)
    savefig(fig, "fig9_monte_carlo.png")


if __name__ == "__main__":
    schematic()
    tradeoff()
    waterfall()
    profiles()
    breakeven()
    tornado()
    monte_carlo()
    print("figures written")
