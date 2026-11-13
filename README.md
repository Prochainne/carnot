# Osmotic energy recovery from desalination brine in the Gulf

Carnot 2026 submission. Model, data, tests and scripts behind the report
"How much of a Gulf desalination plant's electricity can osmotic energy recovery give back?"

## Reproduce everything

```
pip install -r requirements.txt     # tested on Python 3.13; needs numpy >= 2.0;, scipy, pandas, matplotlib, pytest
python run_all.py                   # about 20 minutes on two CPU cores
python run_all.py --quick           # about 10 minutes, 60 Monte Carlo samples instead of 400
```

`run_all.py` runs the tests, then each script in `scripts/` in order. All numbers in the report come
from the files it writes to `results/`; all figures are in `figures/`. The Monte Carlo seed is fixed
(`SEED = 2026` in `scripts/common.py`), so results repeat exactly.

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
| `tests/` | 22 automated checks (`python -m pytest tests`) |
| `data/` | Published data used, with `SOURCES.md` |
| `results/`, `figures/` | Outputs of the last full run |

## Units

SI inside the code (mol/m3, Pa, m/s, m3/s, K). Reports use g/kg, bar, L/m2/h and kWh/m3.
