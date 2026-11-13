"""Checks of the cost arithmetic in scripts/s06_economics.py."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from s06_economics import crf, lcoe, breakeven_membrane_price  # noqa: E402

COST = {"membrane_usd_m2": 10.0, "membrane_share": 0.5, "pretreat_usd_m3d": 400.0,
        "membrane_life_y": 5.0, "discount_rate": 0.07, "om_share": 0.03, "availability": 0.9}


def test_capital_recovery_factor():
    assert crf(0.07, 20) == pytest.approx(0.094393, rel=1e-5)     # textbook annuity value
    assert crf(1e-9, 20) == pytest.approx(1 / 20, rel=1e-6)


def test_lcoe_by_hand():
    # 1 m2 and 1 m3 effluent per m3/day of product, 0.2 kWh/m3 net
    out = lcoe(0.2, 1.0, 1.0, COST)
    k = crf(0.07) + 0.03
    energy = 0.2 * 365 * 0.9
    assert out["lcoe"] == pytest.approx((20 * k + 2 + 400 * k) / energy)
    no_pre = lcoe(0.2, 1.0, 1.0, COST, with_pretreatment=False)
    assert no_pre["lcoe"] == pytest.approx((20 * k + 2) / energy)


def test_no_cost_when_no_net_energy():
    assert np.isnan(lcoe(-0.01, 1.0, 1.0, COST)["lcoe"])


def test_breakeven_price_round_trip():
    price = breakeven_membrane_price(0.2, 3.0, COST, 0.05)
    c = {**COST, "membrane_usd_m2": price}
    assert lcoe(0.2, 3.0, 0.0, c, with_pretreatment=False)["lcoe"] == pytest.approx(0.05)
