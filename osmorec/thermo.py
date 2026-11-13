"""Solution thermodynamics: osmotic pressure and the reversible (ceiling) energy.

All quantities are SI unless stated: concentration c in mol NaCl per m^3 of
solution, pressure in Pa, temperature in K, flows in m^3/s, energy in J.

Salt is represented as NaCl. Sea salt is handled with a single multiplicative
factor on the osmotic pressure (see `SEA_SALT_FACTOR`).
"""
from __future__ import annotations

import numpy as np

R = 8.314462618          # J mol^-1 K^-1
M_NACL = 0.058443        # kg mol^-1
RHO_W = 997.0            # kg m^-3, pure water at 25 C
KWH = 3.6e6              # J per kWh

# Density of NaCl solutions: rho = RHO_W + RHO_SLOPE * c  (c in mol m^-3).
# Slope is the least-squares fit (fixed intercept 997) to the six densities in
# data/nacl_properties_achilli2009.csv; tests/test_thermo.py re-derives it.
RHO_SLOPE = 0.03936      # kg m^-3 per (mol m^-3)

# Pitzer ion-interaction parameters for NaCl(aq) at 25 C (Pitzer & Mayorga 1973).
PITZER_BETA0 = 0.0765
PITZER_BETA1 = 0.2664
PITZER_CPHI = 0.00127
PITZER_APHI = 0.3915     # Debye-Hueckel slope for the osmotic coefficient, 25 C
PITZER_B = 1.2
PITZER_ALPHA = 2.0

# Sea salt carries fewer ions per gram than NaCl: mean molar mass per mole of
# sea-salt particles is 31.40 g/mol (Millero et al. 2008 reference composition)
# against 58.443/2 = 29.22 g/mol for NaCl. At equal mass salinity the osmotic
# pressure is therefore scaled by 29.22/31.40 = 0.93.
SEA_SALT_FACTOR = (M_NACL * 1000 / 2) / 31.40


def density(c):
    """Solution density (kg/m^3) for NaCl concentration c (mol/m^3)."""
    return RHO_W + RHO_SLOPE * np.asarray(c, dtype=float)


def molality(c):
    """Molality (mol per kg water) from molar concentration c (mol/m^3)."""
    c = np.maximum(np.asarray(c, dtype=float), 0.0)
    return c / (density(c) - c * M_NACL)


def conc_from_salinity(s_g_per_kg):
    """Molar concentration (mol/m^3) of an NaCl solution with mass salinity s (g/kg)."""
    w = np.asarray(s_g_per_kg, dtype=float) / 1000.0
    return w * RHO_W / (M_NACL - w * RHO_SLOPE)


def salinity_from_conc(c):
    """Mass salinity (g salt per kg solution) from molar concentration (mol/m^3)."""
    c = np.asarray(c, dtype=float)
    return 1000.0 * c * M_NACL / density(c)


def osmotic_coefficient(m):
    """Pitzer osmotic coefficient of NaCl(aq) at 25 C for molality m (mol/kg)."""
    m = np.maximum(np.asarray(m, dtype=float), 0.0)
    sq = np.sqrt(m)                       # ionic strength I = m for a 1:1 salt
    f_phi = -PITZER_APHI * sq / (1.0 + PITZER_B * sq)
    b_phi = PITZER_BETA0 + PITZER_BETA1 * np.exp(-PITZER_ALPHA * sq)
    return 1.0 + f_phi + m * b_phi + m * m * PITZER_CPHI


def osmotic_pressure(c, T=298.15, ideal=False, factor=1.0):
    """Osmotic pressure (Pa) of a solution with NaCl concentration c (mol/m^3).

    ideal=True  -> van 't Hoff, pi = 2 c R T
    ideal=False -> pi = phi(m) * 2 * m * rho_w * R * T, phi from Pitzer
    `factor` multiplies the result (use SEA_SALT_FACTOR for sea salt).
    """
    c = np.maximum(np.asarray(c, dtype=float), 0.0)
    if ideal:
        return factor * 2.0 * c * R * T
    m = molality(c)
    return factor * osmotic_coefficient(m) * 2.0 * m * RHO_W * R * T


def equilibrium_permeate(QD, cD, QF, cF):
    """Permeate flow at which draw and feed reach equal concentration (no salt leak)."""
    return QD * QF * (cD - cF) / (cD * QD + cF * QF)


def reversible_power(QD, cD, QF, cF, T=298.15, ideal=False, factor=1.0, n=4000):
    """Reversible work rate (W) of mixing a draw (QD, cD) with a feed (QF, cF).

    Water is moved from feed to draw through an ideal semi-permeable membrane
    with the hydraulic pressure always equal to the osmotic pressure
    difference, until the two streams have equal concentration:

        W_rev = integral_0^V*  [pi_D(V) - pi_F(V)] dV

    This equals the Gibbs free energy of mixing released per unit time.
    """
    v_star = equilibrium_permeate(QD, cD, QF, cF)
    v = np.linspace(0.0, v_star, n)
    c_d = cD * QD / (QD + v)
    c_f = cF * QF / (QF - v)
    dpi = (osmotic_pressure(c_d, T, ideal, factor)
           - osmotic_pressure(c_f, T, ideal, factor))
    return float(np.trapezoid(dpi, v))


def counter_current_limit_ideal(cD, cF, T=298.15):
    """Lin, Straub & Elimelech (2014), Eq. 21: maximum specific energy (J per m^3
    of total mixed solution) of constant-pressure counter-current PRO with an
    ideal membrane and van 't Hoff solutions; optimum at dP = dpi/2, phi = 0.5."""
    return 2.0 * R * T / 4.0 * (cD - cF) ** 2 / (cD + cF)
