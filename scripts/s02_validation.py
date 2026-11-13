"""Validation against published data.

1. Osmotic pressure model vs tabulated NaCl data (Achilli et al. 2009).
2. Thermodynamic limits vs the analytical results of Lin, Straub & Elimelech (2014).
3. Flux model, blind prediction: Achilli et al. (2009) with their own membrane parameters.
4. Flux model, calibrate / test split: Straub, Yip & Elimelech (2014) flux data.
5. Plant scale: reported net output of the Fukuoka osmotic power plant.
"""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from common import (DATA, save_json, savefig, plt, BLUE, ORANGE, AQUA, INK2, MUTED, AXIS,
                    BASE_MACH)
from osmorec import (osmotic_pressure, reversible_power, local_flux, solve_module,
                     Membrane, ACHILLI_CTA, CTA_HF, TFC_LAB, LMH, BAR, KWH, Site, Design,
                     evaluate, optimise_pressure, ceiling, SEA_SALT_FACTOR)
from osmorec.thermo import M_NACL, R, counter_current_limit_ideal

out = {}

# 1 --- osmotic pressure ------------------------------------------------------
tab = pd.read_csv(f"{DATA}/nacl_properties_achilli2009.csv", comment="#")
tab = tab[tab.nacl_g_per_L > 0].copy()
c = tab.nacl_g_per_L.values / M_NACL                    # mol/m^3
tab["pitzer_kPa"] = osmotic_pressure(c) / 1e3
tab["vant_hoff_kPa"] = osmotic_pressure(c, ideal=True) / 1e3
tab["pitzer_err_pct"] = 100 * (tab.pitzer_kPa / tab.osmotic_pressure_kPa - 1)
tab["vant_hoff_err_pct"] = 100 * (tab.vant_hoff_kPa / tab.osmotic_pressure_kPa - 1)
tab.to_csv("results/osmotic_pressure_check.csv", index=False)
hi = tab[tab.nacl_g_per_L >= 20]
out["osmotic_pressure"] = {
    "max_abs_err_pct_pitzer_20to60gL": float(hi.pitzer_err_pct.abs().max()),
    "max_abs_err_pct_vant_hoff_20to60gL": float(hi.vant_hoff_err_pct.abs().max()),
    "max_abs_err_kPa_pitzer_all": float((tab.pitzer_kPa - tab.osmotic_pressure_kPa).abs().max()),
}

fig, ax = plt.subplots(figsize=(6.4, 3.8))
cc = np.linspace(0, 1600, 200)
ax.plot(cc * M_NACL, osmotic_pressure(cc) / BAR, color=BLUE, label="This work (Pitzer)")
ax.plot(cc * M_NACL, osmotic_pressure(cc, ideal=True) / BAR, color=MUTED, ls="--", label="van 't Hoff (ideal)")
ax.plot(tab.nacl_g_per_L, tab.osmotic_pressure_kPa / 100, "o", color=ORANGE, ms=7,
        mec="white", mew=1.5, label="Tabulated data (Achilli et al. 2009)", zorder=5)
ax.set_xlabel("NaCl concentration (g/L)")
ax.set_ylabel("Osmotic pressure (bar)")
ax.set_title("Osmotic pressure of NaCl solutions at 25 °C")
ax.legend(loc="upper left")
savefig(fig, "fig1_osmotic_pressure.png")

# 2 --- thermodynamic limits (Lin et al. 2014: 0.256 and 0.192 kWh/m^3) ---------
cD, cF = 600.0, 15.0
phis = np.linspace(0.3, 0.8, 501)
rev = np.array([reversible_power(1 - p, cD, p, cF, ideal=True) for p in phis]) / KWH
dpi = 2 * R * 298.15 * (cD - cF)
ideal_mem = Membrane(A=1e-11, B=0.0, S=0.0, k=np.inf)
cc_num = solve_module(0.5, cD, 0.5, cF, dpi / 2, 1e8, ideal_mem, ideal=True).permeate * dpi / 2 / KWH
out["thermodynamic_limits"] = {
    "reversible_kWh_m3_this_work": float(rev.max()), "reversible_feed_fraction": float(phis[rev.argmax()]),
    "reversible_kWh_m3_Lin2014": 0.256,
    "counter_current_numerical": float(cc_num),
    "counter_current_analytical": float(counter_current_limit_ideal(cD, cF) / KWH),
    "counter_current_kWh_m3_Lin2014": 0.192,
}

# 3 --- Achilli 2009: prediction with the published parameters, nothing fitted --
ach = pd.read_csv(f"{DATA}/achilli2009_pro.csv", comment="#")
pred = []
for _, r in ach.iterrows():
    Jw, _ = local_flux(r.draw_g_per_L / M_NACL, r.feed_g_per_L / M_NACL, r.dP_kPa * 1e3, ACHILLI_CTA)
    pred.append(Jw[0] * 1e6)
ach["Jw_model_1e-6_m_per_s"] = pred
ach["err_pct"] = 100 * (ach["Jw_model_1e-6_m_per_s"] / ach["Jw_1e-6_m_per_s"] - 1)
ach["power_density_model_W_m2"] = ach["Jw_model_1e-6_m_per_s"] * 1e-6 * ach.dP_kPa * 1e3
ach.to_csv("results/validation_achilli.csv", index=False)
out["achilli"] = ach[["draw_g_per_L", "feed_g_per_L", "Jw_1e-6_m_per_s", "Jw_model_1e-6_m_per_s",
                      "err_pct", "power_density_W_m2", "power_density_model_W_m2", "note"]].to_dict("records")

# 4 --- Straub 2014: calibrate on two draw concentrations, test on the other two
st = pd.read_csv(f"{DATA}/straub2014_pro_flux.csv", comment="#")


def flux_lmh(p, d):
    mem = Membrane.from_common_units(*np.exp(p))
    return np.array([local_flux(cM * 1000, 0.0, dp * BAR, mem)[0][0] / LMH
                     for cM, dp in zip(d.cD_mol_per_L, d.dP_bar)])


def fit(d):
    r = least_squares(lambda p: flux_lmh(p, d) - d.Jw_LMH.values, np.log([2.5, 0.5, 500.0]))
    dof = len(d) - 3
    cov = np.linalg.inv(r.jac.T @ r.jac) * (2 * r.cost / dof)
    return r.x, np.sqrt(np.diag(cov))            # log-parameters and their standard errors


def stats(p, d):
    m = flux_lmh(p, d)
    e = m - d.Jw_LMH.values
    return {"n": int(len(d)), "rmse_LMH": float(np.sqrt(np.mean(e ** 2))),
            "mean_abs_err_LMH": float(np.mean(np.abs(e))),
            "bias_LMH": float(np.mean(e)), "mean_measured_LMH": float(d.Jw_LMH.mean())}


train = st[st.cD_mol_per_L.isin([0.6, 2.0])]
test = st[st.cD_mol_per_L.isin([1.0, 3.0])]
p_tr, se_tr = fit(train)
p_all, se_all = fit(st)
# scatter between the two replicate membrane samples at identical conditions
piv = st.pivot_table(index=["cD_mol_per_L", "dP_bar"], columns="sample", values="Jw_LMH")
rep_rms = float(np.sqrt(np.mean((piv[1] - piv[2]) ** 2)))
names = ["A_LMH_per_bar", "B_LMH", "S_um"]
out["straub"] = {
    "calibration_set": "0.6 M and 2 M draw (24 points)", "test_set": "1 M and 3 M draw (26 points)",
    "fitted_on_calibration": dict(zip(names, np.exp(p_tr).tolist())),
    "calibration_stats": stats(p_tr, train), "test_stats": stats(p_tr, test),
    "test_stats_by_draw": {str(cM): stats(p_tr, test[test.cD_mol_per_L == cM]) for cM in [1.0, 3.0]},
    "fitted_on_all": dict(zip(names, np.exp(p_all).tolist())),
    "fitted_on_all_95pct_factor": dict(zip(names, np.exp(1.96 * se_all).tolist())),
    "all_stats": stats(p_all, st),
    "replicate_rms_difference_LMH": rep_rms,
}

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
cols = {0.6: BLUE, 1.0: ORANGE, 2.0: AQUA, 3.0: "#6250d6"}
mem_tr = Membrane.from_common_units(*np.exp(p_tr))
for cM, col in cols.items():
    d = st[st.cD_mol_per_L == cM]
    pp = np.linspace(2, 50 if cM > 1 else (44 if cM == 1 else 25), 80)
    J = np.array([local_flux(cM * 1000, 0.0, x * BAR, mem_tr)[0][0] for x in pp])
    role = "calibration" if cM in (0.6, 2.0) else "test"
    for ax, y, ym in ((axes[0], d.Jw_LMH, J / LMH), (axes[1], d.W_W_m2, J * pp * BAR)):
        ax.plot(pp, ym, color=col, ls="-" if role == "calibration" else "--")
        ax.plot(d.dP_bar, y, "o" if role == "calibration" else "s", color=col, ms=6, mec="white",
                mew=1.2, label=f"{cM:g} M draw ({role})", zorder=5)
axes[0].set_ylabel("Water flux (L m$^{-2}$ h$^{-1}$)")
axes[1].set_ylabel("Power density (W m$^{-2}$)")
axes[0].set_title("Water flux: model lines, measured points")
axes[1].set_title("Power density")
for ax in axes:
    ax.set_xlabel("Applied hydraulic pressure (bar)")
    ax.set_ylim(bottom=0)
axes[1].legend(loc="upper left", fontsize=8.5)
savefig(fig, "fig2_validation_straub.png")

# 5 --- Fukuoka osmotic power plant: 10,000 m^3/d of ~8 % brine + 9,200 m^3/d
#       of treated sewage, reported net output about 110 kW -----------------------
fuk_site = Site(seawater_salinity=35.0, recovery=1 - 35.0 / 80.0, feed_salinity=1.0)
q_brine = 10000.0 / 86400.0
phi = 9200.0 / 10000.0
rows = []
for mem in (CTA_HF, TFC_LAB):
    for sigma in [10, 25, 50, 75, 100, 150, 200, 300, 600, 1500]:
        d, b = optimise_pressure(fuk_site, mem, BASE_MACH, phi, sigma)
        rows.append({"membrane": mem.name, "sigma_m2_per_m3h": sigma, "dP_bar": d.dP / BAR,
                     "net_kW": b.net_per_m3_brine * KWH * q_brine / 1e3,
                     "gross_kW": b.gross_turbine / fuk_site.brine_per_product * KWH * q_brine / 1e3,
                     "power_density_gross_W_m2": b.power_density_gross,
                     "membrane_area_m2": sigma * q_brine * 3600})
fuk = pd.DataFrame(rows)
fuk.to_csv("results/validation_fukuoka.csv", index=False)
ceil_kw = reversible_power(q_brine, fuk_site.c_brine, q_brine * phi, fuk_site.c_feed,
                           factor=SEA_SALT_FACTOR) / 1e3
out["fukuoka"] = {"reported_net_kW": 110.0, "reversible_ceiling_kW": float(ceil_kw),
                  "reported_fraction_of_ceiling": 110.0 / float(ceil_kw), "model": rows}
for mem in (CTA_HF, TFC_LAB):
    f = fuk[fuk.membrane == mem.name]
    out["fukuoka"][mem.name] = {
        "sigma_matching_110kW": float(np.interp(110.0, f.net_kW, f.sigma_m2_per_m3h)),
        "area_matching_110kW_m2": float(np.interp(110.0, f.net_kW, f.membrane_area_m2)),
        "net_kW_large_area": float(f.net_kW.iloc[-1])}

fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9))
ax = axes[0]
m = ach["Jw_1e-6_m_per_s"].max() * 1.15
ax.plot([0, m], [0, m], color=MUTED, lw=1)
ax.fill_between([0, m], [0, 0.85 * m], [0, 1.15 * m], color=MUTED, alpha=0.12, lw=0)
tabd = ach[ach.note.str.startswith("Table")]
rng = ach[~ach.note.str.startswith("Table")]
ax.plot(tabd["Jw_1e-6_m_per_s"], tabd["Jw_model_1e-6_m_per_s"], "o", color=BLUE, ms=8, mec="white",
        mew=1.5, label="Pure-water feed (tabulated)")
ax.plot(rng["Jw_1e-6_m_per_s"], rng["Jw_model_1e-6_m_per_s"], "s", color=ORANGE, ms=7, mec="white",
        mew=1.5, label="5 g/L feed (range quoted in text)")
ax.text(0.97 * m, 0.80 * m, "±15 %", color=INK2, ha="right", fontsize=9)
ax.set_xlabel("Measured flux (10$^{-6}$ m/s)")
ax.set_ylabel("Predicted flux (10$^{-6}$ m/s)")
ax.set_title("Achilli et al. (2009): nothing fitted")
ax.legend(loc="upper left", fontsize=8.5)
ax = axes[1]
for mem, col in ((CTA_HF, BLUE), (TFC_LAB, AQUA)):
    f = fuk[fuk.membrane == mem.name]
    ax.plot(f.sigma_m2_per_m3h, f.net_kW, "-o", color=col, ms=5, label="Model: " + mem.name)
ax.axhline(110, color=ORANGE, lw=2)
ax.text(9, 118, "Reported net output, 110 kW", color=INK2, ha="left", fontsize=9)
ax.set_xscale("log")
ax.legend(loc="lower right", fontsize=8.5)
ax.set_xlabel("Membrane area per brine flow (m$^2$ per m$^3$/h)")
ax.set_ylabel("Net electrical output (kW)")
ax.set_ylim(-60, 260)
ax.axhline(0, color=AXIS, lw=1)
ax.set_title("Fukuoka plant: model vs reported output")
savefig(fig, "fig3_validation_achilli_fukuoka.png")

save_json("validation.json", out)
print(pd.Series(out["thermodynamic_limits"]))
print(ach[["draw_g_per_L", "feed_g_per_L", "Jw_1e-6_m_per_s", "Jw_model_1e-6_m_per_s", "err_pct"]])
print(out["straub"])
print(fuk.round(2).to_string()); print({k: v for k, v in out["fukuoka"].items() if k != "model"})
print(out["osmotic_pressure"])
