"""Checks of the solution thermodynamics against published numbers."""
import os
import numpy as np
import pandas as pd
import pytest

from osmorec import osmotic_pressure, reversible_power, KWH, conc_from_salinity, salinity_from_conc
from osmorec.thermo import (M_NACL, RHO_SLOPE, RHO_W, R, osmotic_coefficient,
                            counter_current_limit_ideal, SEA_SALT_FACTOR)

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
TAB = pd.read_csv(os.path.join(DATA, "nacl_properties_achilli2009.csv"), comment="#")


def test_density_slope_matches_data():
    x = TAB.nacl_g_per_L.values
    y = TAB.density_g_per_L.values - RHO_W
    slope_per_gL = np.sum(x * y) / np.sum(x * x)              # least squares through (0, 997)
    assert RHO_SLOPE == pytest.approx(slope_per_gL * M_NACL * 1000 / 1000, rel=2e-3)


def test_osmotic_pressure_against_table():
    t = TAB[TAB.nacl_g_per_L >= 20]
    model = osmotic_pressure(t.nacl_g_per_L.values / M_NACL) / 1e3
    assert np.max(np.abs(model / t.osmotic_pressure_kPa.values - 1)) < 0.02


def test_pitzer_limits():
    assert osmotic_coefficient(0.0) == pytest.approx(1.0)
    # Robinson & Stokes: phi(NaCl, 1 mol/kg, 25 C) = 0.936
    assert osmotic_coefficient(1.0) == pytest.approx(0.936, abs=0.003)
    # dilute limit tends to van 't Hoff
    c = 1.0
    assert osmotic_pressure(c) == pytest.approx(osmotic_pressure(c, ideal=True), rel=0.03)


def test_salinity_round_trip():
    for s in (1.0, 42.0, 72.4):
        assert float(salinity_from_conc(conc_from_salinity(s))) == pytest.approx(s, rel=1e-12)


def test_sea_salt_factor():
    assert SEA_SALT_FACTOR == pytest.approx(0.93, abs=0.005)


def test_reversible_energy_matches_lin_2014():
    """Seawater (0.6 M) / river water (0.015 M), ideal solutions: 0.256 kWh per m3 mixed."""
    best = max(reversible_power(1 - p, 600.0, p, 15.0, ideal=True) for p in np.linspace(0.5, 0.7, 201))
    assert best / KWH == pytest.approx(0.256, abs=0.001)


def test_reversible_energy_closed_form():
    """Numerical integral equals the ideal Gibbs free energy of mixing."""
    cD, cF, phi = 600.0, 15.0, 0.6
    cM = (1 - phi) * cD + phi * cF
    g = 2 * R * 298.15 * ((1 - phi) * cD * np.log(cD) + phi * cF * np.log(cF) - cM * np.log(cM))
    assert reversible_power(1 - phi, cD, phi, cF, ideal=True, n=20000) == pytest.approx(g, rel=1e-4)


def test_counter_current_formula_value():
    assert counter_current_limit_ideal(600.0, 15.0) / KWH == pytest.approx(0.192, abs=0.001)
