"""Plant-scale energy balance: a seawater reverse-osmosis (SWRO) plant whose
brine is the draw solution of a PRO unit fed with treated sewage effluent (TSE).

Basis: 1 m^3/s of brine. Results are reported per m^3 of desalinated product.
"""
from __future__ import annotations

from dataclasses import dataclass, replace, asdict

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .membrane import Membrane, BAR
from .module import solve_module
from .thermo import (KWH, SEA_SALT_FACTOR, conc_from_salinity, reversible_power,
                     osmotic_pressure, salinity_from_conc)


@dataclass(frozen=True)
class Site:
    """Desalination plant and water sources."""
    seawater_salinity: float = 42.0     # g/kg, Gulf intake
    recovery: float = 0.42              # SWRO water recovery (product / intake)
    sec_ro: float = 3.0                 # kWh per m^3 product, whole SWRO plant
    feed_salinity: float = 1.0          # g/kg, treated sewage effluent
    T: float = 298.15                   # K
    salt_factor: float = SEA_SALT_FACTOR

    @property
    def brine_salinity(self):
        """g/kg; all salt stays in the brine (complete rejection)."""
        return self.seawater_salinity / (1.0 - self.recovery)

    @property
    def brine_per_product(self):
        """m^3 of brine per m^3 of product water."""
        return (1.0 - self.recovery) / self.recovery

    @property
    def c_brine(self):
        return float(conc_from_salinity(self.brine_salinity))

    @property
    def c_feed(self):
        return float(conc_from_salinity(self.feed_salinity))


@dataclass(frozen=True)
class Machines:
    """Rotating equipment, hydraulic losses and pretreatment."""
    eta_turbine: float = 0.85           # hydraulic-to-electric
    eta_pump: float = 0.85              # electric-to-hydraulic
    eta_px: float = 0.95                # pressure exchanger
    dp_feed: float = 1.5 * BAR          # Pa, feed-side friction loss through the module
    dp_draw: float = 0.8 * BAR          # Pa, draw-side friction loss
    e_pre_feed: float = 0.10            # kWh per m^3 of feed (pretreatment + conveyance)
    e_pre_draw: float = 0.0             # kWh per m^3 of brine (already pretreated by the RO plant)


@dataclass(frozen=True)
class Design:
    """Operating point of the PRO unit."""
    dP: float                           # Pa, applied hydraulic pressure difference
    phi: float                          # feed-to-brine volumetric flow ratio Q_F / Q_D
    sigma: float                        # membrane area per brine flow, m^2 per (m^3/h)


@dataclass
class Balance:
    """All terms in kWh per m^3 of desalinated product unless noted."""
    gross_turbine: float
    loss_px: float
    loss_draw_pump: float
    loss_feed_pump: float
    loss_pretreatment: float
    net: float
    fraction_recovered: float           # net / SEC of the RO plant
    net_per_m3_feed: float              # kWh per m^3 TSE used
    net_per_m3_brine: float             # kWh per m^3 brine
    permeate_per_brine: float           # m^3 water transferred per m^3 brine
    feed_utilisation: float             # fraction of the feed that permeates
    power_density_gross: float          # W/m^2 membrane, dP * Jw averaged
    power_density_net: float            # W/m^2 membrane, net electric
    discharge_salinity: float           # g/kg of the diluted brine leaving the PRO unit
    area_per_product: float             # m^2 membrane per (m^3/day) product
    pinched: bool = False

    def as_dict(self):
        return asdict(self)


def evaluate(site: Site, mem: Membrane, mach: Machines, design: Design,
             n_cells=80, ideal=False, flow="counter") -> Balance:
    """Energy balance of the PRO unit for a given design point."""
    QD = 1.0                                          # m^3/s brine (basis)
    QF = design.phi * QD
    area = design.sigma * QD * 3600.0                 # m^2
    m = mem.at_temperature(site.T)
    res = solve_module(QD, site.c_brine, QF, site.c_feed, design.dP, area, m,
                       T=site.T, ideal=ideal, factor=site.salt_factor,
                       n_cells=n_cells, flow=flow)
    dP = design.dP
    # powers in W for 1 m^3/s of brine
    p_turb = mach.eta_turbine * dP * res.permeate
    p_px = (1.0 - mach.eta_px) * dP * QD / mach.eta_pump
    p_draw = mach.dp_draw * QD / mach.eta_pump
    p_feed = mach.dp_feed * QF / mach.eta_pump
    p_pre = (mach.e_pre_feed * QF + mach.e_pre_draw * QD) * KWH
    p_net = p_turb - p_px - p_draw - p_feed - p_pre
    per_brine = 1.0 / (QD * KWH)                      # W -> kWh per m^3 brine
    per_prod = per_brine * site.brine_per_product     # -> kWh per m^3 product
    net = p_net * per_prod
    return Balance(
        gross_turbine=p_turb * per_prod, loss_px=p_px * per_prod,
        loss_draw_pump=p_draw * per_prod, loss_feed_pump=p_feed * per_prod,
        loss_pretreatment=p_pre * per_prod, net=net,
        fraction_recovered=net / site.sec_ro,
        net_per_m3_feed=p_net / (QF * KWH), net_per_m3_brine=p_net * per_brine,
        permeate_per_brine=res.permeate / QD, feed_utilisation=res.permeate / QF,
        power_density_gross=res.power_density, power_density_net=p_net / area,
        discharge_salinity=float(salinity_from_conc(res.cD_out)),
        area_per_product=area * site.brine_per_product / (QD * 86400.0),
        pinched=res.pinched)


def ceiling(site: Site, phi: float, ideal=False) -> float:
    """Reversible (thermodynamic) energy of mixing, kWh per m^3 product."""
    w = reversible_power(1.0, site.c_brine, phi, site.c_feed, site.T, ideal,
                         site.salt_factor)
    return w / KWH * site.brine_per_product


def bulk_osmotic_difference(site: Site) -> float:
    """Osmotic pressure difference between brine and feed at the inlets (Pa)."""
    return float(osmotic_pressure(site.c_brine, site.T, factor=site.salt_factor)
                 - osmotic_pressure(site.c_feed, site.T, factor=site.salt_factor))


def optimise_pressure(site, mem, mach, phi, sigma, n_cells=80):
    """Best applied pressure for fixed feed ratio and membrane area."""
    dpi = bulk_osmotic_difference(site)

    def neg(x):
        return -evaluate(site, mem, mach, Design(x * dpi, phi, sigma), n_cells).net

    r = minimize_scalar(neg, bounds=(0.15, 0.9), method="bounded",
                        options={"xatol": 2e-3})
    d = Design(float(r.x) * dpi, phi, sigma)
    return d, evaluate(site, mem, mach, d, n_cells)


def optimise(site, mem, mach, sigma, n_cells=80, phi_max=3.0, x0=(0.45, 0.9)):
    """Best applied pressure and feed ratio for a fixed membrane area."""
    dpi = bulk_osmotic_difference(site)

    def neg(x):
        f, phi = x
        if not (0.1 <= f <= 0.95 and 0.02 <= phi <= phi_max):
            return 1e3
        return -evaluate(site, mem, mach, Design(f * dpi, phi, sigma), n_cells).net

    r = minimize(neg, x0=np.array(x0), method="Nelder-Mead",
                 options={"xatol": 3e-3, "fatol": 1e-5, "maxiter": 200})
    d = Design(float(r.x[0]) * dpi, float(r.x[1]), sigma)
    return d, evaluate(site, mem, mach, d, n_cells)


def with_(obj, **kw):
    return replace(obj, **kw)


def optimise_feed_value(site, mem, mach, sigma_feed, n_cells=80, x0=(0.5, 0.5)):
    """Best use of a scarce feed: maximise net energy per m^3 of feed over the
    applied pressure and the feed-to-brine ratio, for a fixed membrane area per
    unit of FEED flow, sigma_feed in m^2 per (m^3/h of feed)."""
    dpi = bulk_osmotic_difference(site)

    def neg(x):
        f, phi = x
        if not (0.1 <= f <= 0.95 and 0.02 <= phi <= 3.0):
            return 1e3
        d = Design(f * dpi, phi, sigma_feed * phi)
        return -evaluate(site, mem, mach, d, n_cells).net_per_m3_feed

    r = minimize(neg, x0=np.array(x0), method="Nelder-Mead",
                 options={"xatol": 3e-3, "fatol": 1e-5, "maxiter": 200})
    d = Design(float(r.x[0]) * dpi, float(r.x[1]), sigma_feed * float(r.x[1]))
    return d, evaluate(site, mem, mach, d, n_cells)
