"""Independent cross-check of `module.solve_module`.

Same physical model, different numerics: the membrane area is the independent
variable, split into N cells, and the counter-current boundary-value problem
is solved by under-relaxed fixed-point sweeps on the flux profile. Used only
in the tests; it is slower and does not converge for very large areas.
"""
from __future__ import annotations

import numpy as np

from .membrane import Membrane, local_flux


def solve_module_picard(QD0, cD0, QF0, cF0, dP, area, mem: Membrane, T=298.15,
                        ideal=False, factor=1.0, n_cells=200, relax=0.5,
                        tol=1e-9, max_sweeps=2000):
    """Return (permeate m^3/s, salt leak mol/s, converged)."""
    n = n_cells
    dA = area / n
    Jw = np.zeros(n)
    Js = np.zeros(n)
    converged = False
    for _ in range(max_sweeps):
        cum_w = np.concatenate(([0.0], np.cumsum(Jw * dA)))
        cum_s = np.concatenate(([0.0], np.cumsum(Js * dA)))
        QD = QD0 + cum_w
        nD = QD0 * cD0 - cum_s
        QF = np.maximum(QF0 - (cum_w[-1] - cum_w), 1e-6 * QF0)   # feed enters at node n
        nF = QF0 * cF0 + (cum_s[-1] - cum_s)
        cD = nD / QD
        cF = nF / QF
        Jw_new, Js_new = local_flux(0.5 * (cD[:-1] + cD[1:]), 0.5 * (cF[:-1] + cF[1:]),
                                    dP, mem, T, ideal, factor)
        err = float(np.max(np.abs(Jw_new - Jw))) / max(float(Jw_new.max()), 1e-15)
        Jw = Jw + relax * (Jw_new - Jw)
        Js = Js + relax * (Js_new - Js)
        if err < tol:
            converged = True
            break
    v = float(np.sum(Jw) * dA)
    ok = converged and np.isfinite(v) and 0.0 < v < QF0
    return v, float(np.sum(Js) * dA), bool(ok)
