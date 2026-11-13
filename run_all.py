"""Regenerate every number and figure in the report.

    python run_all.py           full run, about 17 minutes on 16 CPU cores (35-40 on 2)
    python run_all.py --quick   same pipeline with 60 Monte Carlo samples

Outputs go to results/ and figures/. The random seed is fixed in scripts/common.py.
"""
import subprocess
import sys
import time

QUICK = "--quick" in sys.argv
N_MC = "60" if QUICK else "400"
STEPS = [
    ["-m", "pytest", "-q", "tests"],
    ["scripts/s02_validation.py"],
    ["scripts/s03_design.py"],
    ["scripts/s04_uncertainty.py", "cta", "compact", N_MC],
    ["scripts/s04_uncertainty.py", "tfc", "high_recovery", N_MC],
    ["scripts/s04_uncertainty.py", "map"],
    ["scripts/s05_figures.py"],
    # costs and fouling (report sections 8 and 9): they only add files, earlier outputs are unchanged
    ["scripts/s06_economics.py"],
    ["scripts/s07_fouling.py"],
]
t0 = time.time()
for step in STEPS:
    print(">>>", " ".join(step), flush=True)
    subprocess.run([sys.executable, *step], check=True)
print("all steps finished in %.0f min" % ((time.time() - t0) / 60))
