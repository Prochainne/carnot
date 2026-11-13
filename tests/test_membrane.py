"""Checks of the local flux model."""
import numpy as np
import pytest

from osmorec import Membrane, local_flux, ideal_flux, osmotic_pressure, ACHILLI_CTA, CTA_HF, BAR
from osmorec.thermo import M_NACL


def test_reduces_to_ideal_flux_without_polarisation():
    mem = Membrane(A=1e-12, B=0.0, S=0.0, k=np.inf)
    cD, cF, dP = np.array([600.0, 1200.0]), np.array([15.0, 100.0]), 10 * BAR
    Jw, Js = local_flux(cD, cF, dP, mem)
    assert Jw == pytest.approx(ideal_flux(cD, cF, dP, mem), rel=1e-8)
    assert np.all(Js == 0.0)


def test_solution_satisfies_governing_equations():
    mem, cD, cF, dP = CTA_HF, 1200.0, 30.0, 25 * BAR
    Jw, Js = (x[0] for x in local_flux(cD, cF, dP, mem))
    e1, e2 = np.exp(-Jw / mem.k), np.exp(Jw * mem.S / mem.D)
    cDm = cD * e1 - Js / Jw * (1 - e1)
    cFm = cF * e2 + Js / Jw * (e2 - 1)
    assert Js == pytest.approx(mem.B * (cDm - cFm), rel=1e-8)
    assert Jw == pytest.approx(mem.A * (osmotic_pressure(cDm) - osmotic_pressure(cFm) - dP), rel=1e-7)


def test_polarisation_and_leak_reduce_flux():
    base = Membrane(A=1e-12, B=0.0, S=0.0, k=np.inf)
    j0 = local_flux(1000.0, 20.0, 15 * BAR, base)[0][0]
    for kw in ({"S": 5e-4}, {"k": 3e-5}, {"B": 1e-7, "S": 5e-4}):
        m = Membrane(**{**dict(A=1e-12, B=0.0, S=0.0, k=np.inf), **kw})
        assert local_flux(1000.0, 20.0, 15 * BAR, m)[0][0] < j0


def test_flux_decreases_with_pressure_and_vanishes_at_flux_reversal():
    p = np.array([5, 15, 25, 35]) * BAR
    J = np.array([local_flux(800.0, 10.0, x, CTA_HF)[0][0] for x in p])
    assert np.all(np.diff(J) < 0)
    dpi = float(osmotic_pressure(800.0) - osmotic_pressure(10.0))
    assert local_flux(800.0, 10.0, 1.01 * dpi, CTA_HF)[0][0] == 0.0


def test_achilli_2009_blind_prediction():
    """Tabulated fluxes at 970 kPa with pure-water feed, published A, B, S, k; nothing fitted."""
    for draw_gL, measured in ((35.0, 2.81e-6), (60.0, 5.21e-6)):
        Jw = local_flux(draw_gL / M_NACL, 0.0, 970e3, ACHILLI_CTA)[0][0]
        assert Jw == pytest.approx(measured, rel=0.05)


def test_temperature_scaling_direction():
    warm = CTA_HF.at_temperature(308.15)
    assert warm.A > CTA_HF.A and warm.D > CTA_HF.D and warm.B > CTA_HF.B
    assert CTA_HF.at_temperature(298.15) is CTA_HF
