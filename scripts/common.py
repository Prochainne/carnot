"""Shared settings: base-case parameters, uncertainty ranges, plotting style."""
import json
import os
import sys
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from osmorec import Site, Machines, BAR  # noqa: E402

DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(FIGURES, exist_ok=True)

SEED = 2026
BASE_SITE = Site()
BASE_MACH = Machines()

# name: (low, base, high, label). Membrane entries are multipliers on the base membrane.
RANGES = {
    "seawater_salinity": (40.0, 42.0, 45.0, "Seawater salinity (g/kg)"),
    "recovery":          (0.38, 0.42, 0.48, "RO recovery"),
    "sec_ro":            (2.8, 3.0, 4.0, "RO energy use (kWh/m3)"),
    "feed_salinity":     (0.5, 1.0, 2.5, "Effluent salinity (g/kg)"),
    "T":                 (293.15, 298.15, 308.15, "Temperature (20-35 C)"),
    "salt_factor":       (0.90, 0.93, 0.96, "Sea-salt osmotic factor"),
    "A_mult":            (0.5, 1.0, 2.0, "Water permeability A (x0.5-x2)"),
    "B_mult":            (0.5, 1.0, 2.0, "Salt permeability B (x0.5-x2)"),
    "S_mult":            (0.5, 1.0, 1.5, "Structural parameter S (x0.5-x1.5)"),
    "k_mult":            (0.5, 1.0, 2.0, "Draw mass transfer k (x0.5-x2)"),
    "eta_turbine":       (0.80, 0.85, 0.90, "Turbine efficiency"),
    "eta_pump":          (0.75, 0.85, 0.90, "Pump efficiency"),
    "eta_px":            (0.92, 0.95, 0.97, "Pressure-exchanger efficiency"),
    "dp_feed":           (0.5 * BAR, 1.5 * BAR, 4.0 * BAR, "Feed-side pressure loss (0.5-4 bar)"),
    "dp_draw":           (0.4 * BAR, 0.8 * BAR, 1.5 * BAR, "Draw-side pressure loss (0.4-1.5 bar)"),
    "e_pre_feed":        (0.05, 0.10, 0.30, "Effluent pretreatment (kWh/m3)"),
    "e_pre_draw":        (0.0, 0.0, 0.05, "Brine pretreatment (kWh/m3)"),
}
SITE_KEYS = {"seawater_salinity", "recovery", "sec_ro", "feed_salinity", "T", "salt_factor"}
MEM_KEYS = {"A_mult": "A", "B_mult": "B", "S_mult": "S", "k_mult": "k"}


def build(params, mem_base, site=BASE_SITE, mach=BASE_MACH):
    """Apply a dict of parameter values to the base site, membrane and machines."""
    s_kw = {k: v for k, v in params.items() if k in SITE_KEYS}
    m_kw = {MEM_KEYS[k]: getattr(mem_base, MEM_KEYS[k]) * v for k, v in params.items() if k in MEM_KEYS}
    x_kw = {k: v for k, v in params.items() if k not in SITE_KEYS and k not in MEM_KEYS}
    return replace(site, **s_kw), replace(mem_base, **m_kw), replace(mach, **x_kw)


def save_json(name, obj):
    with open(os.path.join(RESULTS, name), "w") as f:
        json.dump(obj, f, indent=2, default=float)


# --- plotting style (one axis per chart, thin marks, recessive grid) -----------
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 10, "text.color": INK,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "lines.linewidth": 2.0, "legend.frameon": False, "savefig.dpi": 200, "savefig.bbox": "tight",
})


def savefig(fig, name):
    fig.savefig(os.path.join(FIGURES, name))
    plt.close(fig)
