# Osmotic energy recovery from desalination brine in the Gulf

[![tests](https://github.com/Prochainne/carnot/actions/workflows/tests.yml/badge.svg)](https://github.com/Prochainne/carnot/actions/workflows/tests.yml)

Carnot 2026 submission by team Finners. Model, data, tests and scripts behind the report
"How much of a Gulf desalination plant's electricity can osmotic energy recovery give back?"
([PDF](report/gulf_osmotic_energy_recovery_report.pdf)).

## Reproduce everything

```
pip install -r requirements.txt     # tested on Python 3.13 and 3.14; needs numpy >= 2.0, scipy, pandas, matplotlib, pytest
python run_all.py                   # about 17 minutes on 16 CPU cores, 35-40 minutes on 2
python run_all.py --quick           # 60 Monte Carlo samples instead of 400
```

The scripts use every CPU core available. The longest steps are the design study (`s03_design.py`,
about 6 minutes, which has only two membranes to work on in parallel) and the fouling Monte Carlo
(`s07_fouling.py`, about 5 minutes).

`run_all.py` runs the tests, then each script in `scripts/` in order. All numbers in the report come
from the files it writes to `results/`; figure N of the report is `figures/figN_*.png`. The Monte Carlo
seed is fixed (`SEED = 2026` in `scripts/common.py`), so results repeat exactly.

## Key results

Share of the desalination plant's electricity recovered, base case, from `results/designs.json`
(commercial CTA hollow-fibre membrane) and `results/uncertainty_*.json`:

| Case | Recovered fraction |
|---|---|
| Compact design (CTA) | 6.0 % |
| High-recovery design (CTA) | 15.8 % |
| Large-area design (CTA) | 19.7 % |
| Monte Carlo, CTA compact, 400 samples: median (5th-95th percentile) | 3.9 % (0.9-7.3 %) |
| Monte Carlo, TFC high-recovery, 400 samples: median (5th-95th percentile) | 6.3 % (1.3-11.5 %) |

Net energy is positive in about 98 % of Monte Carlo samples. Effluent pretreatment energy and membrane
water permeability drive most of the uncertainty.

In short: about 4 to 6 % of the plant's electricity comes back in realistic conditions, at a cost of
about US$1 per kWh, some 70 times the price of Gulf solar power. In Abu Dhabi the effluent supply, not
the technology, is the binding limit: with today's unused effluent the Taweelah plant would recover
1 to 3 %.

## Costs and fouling (report sections 8 and 9)

These two analyses were added after the first draft of the report. They only add files: the results
of sections 6 and 7 are unchanged. Both are first estimates, with assumptions listed in
`data/SOURCES.md`.

**Costs (`s06_economics.py`, Figures 10 and 11).** Cost of each kWh the PRO unit makes, from membrane
price and life, the rest of the PRO equipment, effluent pretreatment plant, discount rate, O&M and
availability.

| Case (base-case costs) | US$ per kWh | Without pretreatment capital |
|---|---|---|
| CTA compact | 0.98 | 0.44 |
| CTA high-recovery | 1.24 | 0.78 |
| TFC high-recovery | 0.95 | 0.26 |
| Monte Carlo, CTA compact: median (5th-95th percentile) | 1.9 (0.9-7.3) | 0.9 (0.3-4.5) |

For comparison, Abu Dhabi's Al Dhafra solar plant bid US$0.0135 per kWh. No design or Monte Carlo
sample comes close: the membrane would have to be nearly free (under US$1/m2) even without any
pretreatment cost. The CTA compact cost without pretreatment, US$0.44/kWh, matches the minimum
Chung et al. (2017) found for RO brine and wastewater. Recovering this energy is technically possible
but far from paying for itself at today's membrane and pretreatment costs.

The cost tornado (Figure 11) shows which cost input matters most: effluent pretreatment capital and
membrane price, well ahead of discount rate, membrane life, O&M, availability and equipment share.
Even the cheapest end of any single input leaves the cost above US$0.50/kWh.

**Fouling (`s07_fouling.py`, Figure 12).** The main results assume clean membranes. Treated sewage can coat the
membrane and clog its support layer; with the time-averaged fouling severity set from 0 to 40 %, the CTA
compact design falls from 6.0 % to 3.4-3.6 % of plant electricity and the TFC high-recovery design from
11.2 % to 8.0-9.7 %. The fouling levels are illustrative: no PRO pilot data on Gulf effluent was found.

**Everything together (`s07_fouling.py`, `results/fouling_monte_carlo.json`).** The report's 400 Monte
Carlo samples rerun with fouling added (severity 0-40 %, most likely 15 %, coating + support) and
costs sampled as above:

| Design | Plant electricity recovered: median (5th-95th) | US$ per kWh: median (5th-95th) |
|---|---|---|
| CTA compact | 2.9 % (0.1-6.0 %) | 2.4 (1.0-72) |
| TFC high-recovery | 5.3 % (0.3-10.0 %) | 2.2 (1.0-28) |

The very high upper costs are samples where the unit barely makes any net power. Net power is positive
in 96 % of samples; no sample costs under US$0.10/kWh.

## Layout

| Path | What it is |
|---|---|
| `osmorec/thermo.py` | Osmotic pressure (Pitzer), density, reversible energy of mixing |
| `osmorec/membrane.py` | Local water and salt flux with concentration polarisation |
| `osmorec/module.py` | Counter-current membrane stage of finite area |
| `osmorec/module_picard.py` | Second, independent numerical scheme, used only as a cross-check |
| `osmorec/system.py` | Plant energy balance and optimisation of the operating point |
| `scripts/common.py` | Base-case inputs, uncertainty ranges, plot style |
| `scripts/s02_validation.py` | Comparison with published data (Figures 1-3) |
| `scripts/s03_design.py` | Optimum design, waterfall, break-even, Abu Dhabi scenarios |
| `scripts/s04_uncertainty.py` | Sensitivity, Monte Carlo, break-even map |
| `scripts/s05_figures.py` | Draws Figures 0 and 4-9 from `results/` |
| `scripts/s06_economics.py` | Cost per kWh, cost tornado, cost Monte Carlo (Figures 10 and 11) |
| `scripts/s07_fouling.py` | Fouled membranes, and a Monte Carlo with fouling and costs together (Figure 12) |
| `tests/` | 26 automated checks (`python -m pytest tests`), also run on GitHub for every push |
| `data/` | Published data used, with `SOURCES.md` |
| `results/`, `figures/` | Outputs of the last full run |
| `report/` | The submitted report (PDF) |

## Units

SI inside the code (mol/m3, Pa, m/s, m3/s, K). Reports use g/kg, bar, L/m2/h and kWh/m3.
