"""Checks of the module solver and the plant energy balance."""
from dataclasses import replace

import numpy as np
import pytest

from osmorec import (Membrane, solve_module, CTA_HF, TFC_LAB, Site, Machines, Design, evaluate,
                     ceiling, bulk_osmotic_difference, KWH, BAR)
from osmorec.module_picard import solve_module_picard
from osmorec.thermo import R, counter_current_limit_ideal

SITE, MACH = Site(), Machines()
ARGS = dict(T=SITE.T, factor=SITE.salt_factor)


def test_ideal_large_area_limit_matches_analytical_result():
    """Lin et al. (2014) Eq. 21: ideal membrane, van 't Hoff, dP = dpi/2, equal flows."""
    cD, cF = 600.0, 15.0
    dpi = 2 * R * 298.15 * (cD - cF)
    ideal = Membrane(A=1e-11, B=0.0, S=0.0, k=np.inf)
    r = solve_module(0.5, cD, 0.5, cF, dpi / 2, 1e8, ideal, ideal=True)
    assert r.pinched
    assert r.permeate * dpi / 2 == pytest.approx(counter_current_limit_ideal(cD, cF), rel=1e-5)


def test_two_independent_numerical_schemes_agree():
    for sigma in (50.0, 150.0):
        a = solve_module(1.0, SITE.c_brine, 0.8, SITE.c_feed, 28 * BAR, sigma * 3600, CTA_HF,
                         n_cells=400, **ARGS)
        v, leak, ok = solve_module_picard(1.0, SITE.c_brine, 0.8, SITE.c_feed, 28 * BAR,
                                          sigma * 3600, CTA_HF, **ARGS)
        assert ok
        assert a.permeate == pytest.approx(v, rel=1e-4)
        assert a.salt_leak == pytest.approx(leak, rel=1e-3)


def test_grid_convergence():
    coarse = solve_module(1.0, SITE.c_brine, 0.8, SITE.c_feed, 28 * BAR, 150 * 3600, CTA_HF, n_cells=80, **ARGS)
    fine = solve_module(1.0, SITE.c_brine, 0.8, SITE.c_feed, 28 * BAR, 150 * 3600, CTA_HF, n_cells=640, **ARGS)
    assert coarse.permeate == pytest.approx(fine.permeate, rel=1e-4)


def test_mass_balances_close():
    r = solve_module(1.0, SITE.c_brine, 0.8, SITE.c_feed, 28 * BAR, 150 * 3600, TFC_LAB, n_cells=400, **ARGS)
    assert r.water_balance_error < 1e-9
    assert r.salt_balance_error < 1e-4
    salt_in = 1.0 * SITE.c_brine + 0.8 * SITE.c_feed
    salt_out = r.QD_out * r.cD_out + r.QF_out * r.cF_out
    assert salt_out == pytest.approx(salt_in, rel=1e-12)
    assert r.QD_out + r.QF_out == pytest.approx(1.8, rel=1e-12)


def test_more_area_never_reduces_permeate_and_counter_current_beats_co_current():
    v = [solve_module(1.0, SITE.c_brine, 0.8, SITE.c_feed, 28 * BAR, s * 3600, CTA_HF, **ARGS).permeate
         for s in (20, 60, 180, 540)]
    assert np.all(np.diff(v) > 0)
    co = solve_module(1.0, SITE.c_brine, 0.8, SITE.c_feed, 28 * BAR, 540 * 3600, CTA_HF, flow="co", **ARGS)
    assert co.permeate < v[-1]


def test_energy_balance_adds_up_and_respects_the_ceiling():
    d = Design(26 * BAR, 0.9, 200.0)
    b = evaluate(SITE, CTA_HF, MACH, d)
    assert b.net == pytest.approx(b.gross_turbine - b.loss_px - b.loss_draw_pump
                                  - b.loss_feed_pump - b.loss_pretreatment, rel=1e-12)
    assert b.fraction_recovered == pytest.approx(b.net / SITE.sec_ro)
    assert b.gross_turbine / MACH.eta_turbine < ceiling(SITE, d.phi)      # second law
    assert SITE.seawater_salinity < b.discharge_salinity < SITE.brine_salinity


def test_losses_reduce_net_energy_monotonically():
    d = Design(26 * BAR, 0.9, 200.0)
    base = evaluate(SITE, CTA_HF, MACH, d).net
    for kw in ({"e_pre_feed": 0.3}, {"dp_feed": 4 * BAR}, {"eta_turbine": 0.7}, {"eta_px": 0.9}):
        assert evaluate(SITE, CTA_HF, replace(MACH, **kw), d).net < base


def test_site_bookkeeping():
    assert SITE.brine_salinity == pytest.approx(42.0 / 0.58)
    assert SITE.brine_per_product == pytest.approx(0.58 / 0.42)
    assert 50 * BAR < bulk_osmotic_difference(SITE) < 65 * BAR
