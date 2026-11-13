"""Module-scale model: a membrane stage of finite area with the draw and feed
streams in plug flow, counter-current (or co-current), at constant applied
pressure.

Steady-state balances along the membrane area a:

    dQ_D/da = +Jw     d(Q_D c_D)/da = -Js      (draw gains water, loses salt)
    dQ_F/da = -Jw     d(Q_F c_F)/da = +Js      (feed, along its own direction)

Counter-current flow is a two-point boundary-value problem in a. We avoid
solving it as such by using the draw flow q = Q_D as the independent variable.
For a trial total permeate V, overall balances fix every bulk state as a
function of q alone (q runs from Q_D0 to Q_D0 + V), so the local flux Jw(q)
can be evaluated for all positions at once, and the membrane area needed to
transfer V is

    a_req(V) = integral_{Q_D0}^{Q_D0+V} dq / Jw(q).

a_req increases monotonically with V and diverges when the driving force
vanishes somewhere (the thermodynamic pinch). The permeate of a stage with
area `area` is the root of a_req(V) = area, found with Brent's method.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .membrane import Membrane, local_flux


@dataclass
class ModuleResult:
    permeate: float          # m^3/s of water moved from feed to draw
    salt_leak: float         # mol/s of salt moved from draw to feed
    cD_out: float            # mol/m^3
    cF_out: float            # mol/m^3
    QD_out: float            # m^3/s
    QF_out: float            # m^3/s
    power_density: float     # W/m^2, dP * permeate / area (gross, membrane level)
    area_profile: np.ndarray  # cumulative area at each cell edge, from the draw inlet
    Jw: np.ndarray           # m/s, per cell (cell 0 = draw inlet)
    cD: np.ndarray           # bulk concentrations at cell centres
    cF: np.ndarray
    pinched: bool            # True if the stage is limited by a vanishing driving force
    water_balance_error: float
    salt_balance_error: float


def _profiles(V, QD0, cD0, QF0, cF0, dP, mem, T, ideal, factor, n, flow, passes=3):
    """Bulk states and fluxes on a midpoint grid in q for a trial permeate V."""
    dq = V / n
    q = QD0 + (np.arange(n) + 0.5) * dq            # draw flow at cell centres
    w = q - QD0                                    # water gained by the draw so far
    leak_mid = np.zeros(n)                         # salt lost by the draw up to cell centre
    leak_tot = 0.0
    Jw = Js = None
    for _ in range(passes):
        nD = QD0 * cD0 - leak_mid
        if flow == "counter":
            QF = QF0 - V + w                       # feed still carries what has yet to permeate
            nF = QF0 * cF0 + (leak_tot - leak_mid)
        else:                                      # co-current
            QF = QF0 - w
            nF = QF0 * cF0 + leak_mid
        cD = np.maximum(nD, 0.0) / q
        cF = np.maximum(nF, 0.0) / QF
        Jw, Js = local_flux(cD, cF, dP, mem, T, ideal, factor)
        if mem.B == 0.0 or np.any(Jw <= 0.0):
            break
        ratio = Js / Jw                            # mol of salt per m^3 of permeate
        cum_edges = np.concatenate(([0.0], np.cumsum(ratio * dq)))
        leak_mid = 0.5 * (cum_edges[:-1] + cum_edges[1:])
        leak_tot = cum_edges[-1]
    return q, dq, Jw, Js, cD, cF, leak_tot


def required_area(V, QD0, cD0, QF0, cF0, dP, mem, T=298.15, ideal=False,
                  factor=1.0, n=80, flow="counter"):
    """Membrane area needed to transfer a permeate V (inf if not attainable)."""
    _, dq, Jw, _, _, _, leak = _profiles(V, QD0, cD0, QF0, cF0, dP, mem, T, ideal, factor, n, flow)
    if np.any(Jw <= 0.0):
        return np.inf
    # the driving force must also be positive at the two ends of the stage,
    # which the midpoint grid does not sample
    nD_out = QD0 * cD0 - leak
    if flow == "counter":
        ends_D = np.array([cD0, nD_out / (QD0 + V)])
        ends_F = np.array([(QF0 * cF0 + leak) / (QF0 - V), cF0])
    else:
        ends_D = np.array([cD0, nD_out / (QD0 + V)])
        ends_F = np.array([cF0, (QF0 * cF0 + leak) / (QF0 - V)])
    J_ends, _ = local_flux(ends_D, ends_F, dP, mem, T, ideal, factor)
    if np.any(J_ends <= 0.0):
        return np.inf
    return float(np.sum(dq / Jw))


def solve_module(QD0, cD0, QF0, cF0, dP, area, mem: Membrane, T=298.15,
                 ideal=False, factor=1.0, n_cells=80, flow="counter") -> ModuleResult:
    """Steady state of one PRO stage.

    QD0, cD0 : draw inlet flow (m^3/s) and concentration (mol/m^3)
    QF0, cF0 : feed inlet flow and concentration
    dP       : applied hydraulic pressure difference (Pa), uniform over the stage
    area     : membrane area (m^2)
    """
    args = (QD0, cD0, QF0, cF0, dP, mem, T, ideal, factor, n_cells, flow)

    def h(V):                                      # > 0 while the area is not used up
        return area / required_area(V, *args[:-2], n=n_cells, flow=flow) - 1.0

    v_hi = QF0 * (1.0 - 1e-9)
    v_lo = QF0 * 1e-9
    pinched = False
    if h(v_lo) <= 0.0:                             # no positive flux at the inlet
        V = 0.0
    elif h(v_hi) > 0.0:                            # the whole feed permeates
        V = v_hi
    else:
        V = brentq(h, v_lo, v_hi, xtol=1e-10 * QF0, rtol=1e-10, maxiter=200)
        # the stage is pinch-limited if slightly more permeate is unattainable
        pinched = not np.isfinite(required_area(min(V * 1.001, v_hi), *args[:-2],
                                                n=n_cells, flow=flow))
    if V <= 0.0:
        z = np.zeros(n_cells)
        return ModuleResult(0.0, 0.0, cD0, cF0, QD0, QF0, 0.0, np.linspace(0, area, n_cells + 1),
                            z, z + cD0, z + cF0, False, 0.0, 0.0)
    q, dq, Jw, Js, cD, cF, leak = _profiles(V, *args)
    a_edges = np.concatenate(([0.0], np.cumsum(dq / Jw)))
    QD_out = QD0 + V
    QF_out = QF0 - V
    nD_out = QD0 * cD0 - leak
    nF_out = QF0 * cF0 + leak
    # balance checks: compare flux integrals over the area with the stream changes
    water_from_flux = float(np.sum(Jw * np.diff(a_edges)))
    salt_from_flux = float(np.sum(Js * np.diff(a_edges)))
    water_err = abs(water_from_flux - V) / V
    salt_err = abs(salt_from_flux - leak) / max(QD0 * cD0, 1e-300)
    return ModuleResult(
        permeate=float(V), salt_leak=float(leak),
        cD_out=float(nD_out / QD_out), cF_out=float(nF_out / QF_out),
        QD_out=float(QD_out), QF_out=float(QF_out),
        power_density=dP * V / area,
        area_profile=a_edges, Jw=Jw, cD=cD, cF=cF, pinched=bool(pinched),
        water_balance_error=float(water_err), salt_balance_error=float(salt_err))
