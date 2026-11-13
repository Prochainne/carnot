"""Local (coupon-scale) water and salt flux across a PRO membrane.

Film-theory model with external concentration polarisation on the draw side,
internal concentration polarisation in the porous support on the feed side and
reverse salt flux, following Yip et al. (2011) but keeping the osmotic pressure
non-linear in concentration (so interface concentrations are computed
explicitly instead of interface osmotic pressures).

Orientation: active layer faces the draw solution (PRO mode).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .thermo import osmotic_pressure

LMH = 1.0 / 3.6e6            # 1 L m^-2 h^-1 in m/s
BAR = 1.0e5                  # Pa
D_NACL_25 = 1.49e-9          # m^2/s, NaCl in water at 25 C (Achilli et al. 2009, Table)
T_REF = 298.15
C_SATURATION = 5400.0        # mol/m^3, about the solubility of NaCl in water at 25 C


@dataclass(frozen=True)
class Membrane:
    """Intrinsic membrane and hydrodynamic parameters at 25 C (SI units)."""
    A: float                 # water permeability, m s^-1 Pa^-1
    B: float                 # salt permeability, m s^-1
    S: float                 # support structural parameter, m
    k: float = 38.5e-6       # draw-side mass transfer coefficient, m s^-1
    D: float = D_NACL_25     # salt diffusivity, m^2 s^-1
    name: str = "membrane"

    @classmethod
    def from_common_units(cls, A_lmh_bar, B_lmh, S_um, k_um_s=38.5, D=D_NACL_25, name="membrane"):
        return cls(A=A_lmh_bar * LMH / BAR, B=B_lmh * LMH, S=S_um * 1e-6,
                   k=k_um_s * 1e-6, D=D, name=name)

    def at_temperature(self, T):
        """Scale parameters from 25 C to T (K).

        A ~ 1/mu (viscous flow of water), D ~ T/mu (Stokes-Einstein),
        B ~ D (salt diffuses through the active layer), k ~ D^(2/3)
        (Sherwood correlation at fixed hydrodynamics).
        """
        if abs(T - T_REF) < 1e-9:
            return self
        mu_ratio = water_viscosity(T_REF) / water_viscosity(T)
        d_ratio = (T / T_REF) * mu_ratio
        return replace(self, A=self.A * mu_ratio, B=self.B * d_ratio,
                       D=self.D * d_ratio, k=self.k * d_ratio ** (2.0 / 3.0))


def water_viscosity(T):
    """Dynamic viscosity of water (Pa s), T in K. Standard Vogel-type correlation."""
    return 2.414e-5 * 10.0 ** (247.8 / (T - 140.0))


def local_flux(cD, cF, dP, mem: Membrane, T=298.15, ideal=False, factor=1.0,
               rtol=1e-10, maxit=80):
    """Water flux Jw (m/s) and reverse salt flux Js (mol m^-2 s^-1).

    cD, cF : bulk draw / feed concentrations (mol/m^3), scalars or arrays
    dP     : applied hydraulic pressure difference (Pa)

    Solves   Jw = A [ pi(cD_m) - pi(cF_m) - dP ]   with interface concentrations

        cD_m = cD exp(-Jw/k)     - (Js/Jw) (1 - exp(-Jw/k))
        cF_m = cF exp(Jw S/D)    + (Js/Jw) (exp(Jw S/D) - 1)
        Js   = B (cD_m - cF_m)

    by a bracketed false-position (Illinois) iteration, vectorised over cells.
    If the polarised driving force is not positive, Jw = 0 (no back-flow is
    modelled) and salt still diffuses.
    """
    cD = np.atleast_1d(np.asarray(cD, dtype=float))
    cF = np.atleast_1d(np.asarray(cF, dtype=float))
    cD, cF = np.broadcast_arrays(cD, cF)
    A, B, S, k, D = mem.A, mem.B, mem.S, mem.k, mem.D
    inv_k = 0.0 if not np.isfinite(k) else 1.0 / k
    K = S / D

    def pi(c):
        return osmotic_pressure(c, T, ideal, factor)

    def resid(J):
        e1 = np.exp(-J * inv_k)
        e2 = np.exp(np.minimum(J * K, 50.0))
        dcm = (cD * e1 - cF * e2) / (1.0 + (B / J) * (e2 - e1))
        Js = B * dcm
        cDm = cD * e1 - (Js / J) * (1.0 - e1)
        # the support cannot hold more salt than a saturated NaCl solution
        cFm = np.minimum(cF * e2 + (Js / J) * (e2 - 1.0), C_SATURATION)
        return A * (pi(cDm) - pi(cFm) - dP) - J, Js

    tiny = 1e-13
    j_max = A * (pi(cD) - pi(cF) - dP)
    lo = np.full(cD.shape, tiny)
    f_lo, js_lo = resid(lo)
    active = (j_max > 2 * tiny) & (f_lo > 0.0)
    hi = np.where(active, j_max, 2 * tiny)
    f_hi, _ = resid(hi)
    f_hi = np.where(active, np.minimum(f_hi, -1e-30), -1.0)
    f_lo = np.where(active, f_lo, 1.0)
    side = np.zeros(cD.shape, dtype=int)
    scale = max(float(np.max(np.abs(j_max))), tiny)
    x = lo.copy()
    for _ in range(maxit):
        x = (lo * f_hi - hi * f_lo) / (f_hi - f_lo)
        fx, _ = resid(x)
        pos = fx > 0.0
        f_hi = np.where(pos & (side == 1), 0.5 * f_hi, f_hi)
        f_lo = np.where(~pos & (side == -1), 0.5 * f_lo, f_lo)
        lo = np.where(pos, x, lo)
        f_lo = np.where(pos, fx, f_lo)
        hi = np.where(~pos, x, hi)
        f_hi = np.where(~pos, fx, f_hi)
        side = np.where(pos, 1, -1)
        if np.max(np.abs(fx[active]), initial=0.0) < rtol * scale:
            break
    Jw = np.where(active, x, 0.0)
    _, js_x = resid(np.where(active, x, tiny))
    Js = np.where(active, js_x, js_lo)
    return Jw, Js


def ideal_flux(cD, cF, dP, mem: Membrane, T=298.15, ideal=False, factor=1.0):
    """Flux with no polarisation and no salt leak: Jw = A (dpi - dP), clipped at 0."""
    dpi = osmotic_pressure(cD, T, ideal, factor) - osmotic_pressure(cF, T, ideal, factor)
    return np.maximum(mem.A * (dpi - dP), 0.0)
