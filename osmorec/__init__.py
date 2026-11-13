"""osmorec: osmotic energy recovery from desalination brine (Carnot 2026)."""
from .thermo import (osmotic_pressure, osmotic_coefficient, conc_from_salinity,
                     salinity_from_conc, reversible_power, SEA_SALT_FACTOR, KWH)
from .membrane import Membrane, local_flux, ideal_flux, LMH, BAR
from .module import solve_module
from .system import (Site, Machines, Design, Balance, evaluate, ceiling,
                     optimise, optimise_pressure, optimise_feed_value,
                     bulk_osmotic_difference)

# Membranes used in the study -------------------------------------------------
# Commercial cellulose-triacetate hollow fibre (Toyobo module 1 in Tanaka et al. 2018)
CTA_HF = Membrane.from_common_units(0.27, 0.035, 1024, name="CTA hollow fibre (commercial)")
# Flat-sheet CTA of Achilli et al. (2009): A, B, S and k as reported there
ACHILLI_CTA = Membrane(A=1.87e-12, B=1.11e-7, S=6.78e-4, k=8.48e-5, D=1.49e-9,
                       name="CTA flat sheet (Achilli 2009)")
# Lab-scale thin-film composite: parameters fitted in this study to all 50 flux
# measurements of Straub, Yip & Elimelech (2014); see scripts/s02_validation.py
TFC_LAB = Membrane.from_common_units(2.72, 0.995, 459, name="TFC flat sheet (lab, fitted)")
