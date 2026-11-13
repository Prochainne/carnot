# Data sources

Every number that enters the model from outside is listed here with where it came from.

## Files in this folder

| File | Content | Source |
|---|---|---|
| `nacl_properties_achilli2009.csv` | Osmotic pressure, viscosity, density, diffusivity of NaCl solutions at 25 C | Achilli, Cath & Childress, J. Membr. Sci. 343 (2009) 42-52, as printed in A. Achilli's PhD dissertation (Univ. Nevada, Reno, 2009), https://scholarwolf.unr.edu/server/api/core/bitstreams/1433357a-f588-42f4-8ba8-695faf2f5016/content |
| `achilli2009_pro.csv` | Measured PRO water flux and power density at 970 kPa | Same source |
| `straub2014_pro_flux.csv` | 50 measured PRO fluxes, thin-film composite membrane, 0.6-3 M NaCl draw | Straub, Yip & Elimelech, Environ. Sci. Technol. Lett. 1 (2014) 55-59, Supporting Information Table S1, https://core.ac.uk/download/161454577.pdf |
| `plants.csv` | Plant-scale figures for Taweelah, Fukuoka and Abu Dhabi effluent | See links below |

## Parameters set in code

| Parameter | Value | Source |
|---|---|---|
| Pitzer parameters for NaCl at 25 C | beta0 0.0765, beta1 0.2664, Cphi 0.00127, Aphi 0.3915 | Pitzer & Mayorga, J. Phys. Chem. 77 (1973) 2300 (standard textbook values) |
| Mean molar mass of sea-salt particles | 31.40 g/mol | Nayar et al., Desalination 390 (2016) 1-24, after Millero et al. (2008) |
| Commercial CTA hollow-fibre membrane | A 0.27 L/m2/h/bar, B 0.035 L/m2/h, S 1024 um | Tanaka et al., Desalination (2018), Toyobo module (1), https://da.lib.kobe-u.ac.jp/da/kernel/90004891/90004891.pdf |
| Flat-sheet CTA membrane (validation) | A 1.87e-9 m/s/kPa, B 1.11e-7 m/s, S 6.78e-4 m, k 8.48e-5 m/s | Achilli et al. (2009) |
| Draw-side mass transfer coefficient | 38.5 um/s | Straub, Deshmukh & Elimelech, Energy Environ. Sci. 9 (2016) 31-48 |
| Thermodynamic limits used as checks | 0.256 and 0.192 kWh/m3 | Lin, Straub & Elimelech, Energy Environ. Sci. 7 (2014) 2706-2714 |
| Pretreatment energy range, component efficiencies, pressure losses | 0.05-0.2 kWh/m3; turbine up to 90 %; 0.8 bar/m draw, 2-5 bar/m feed | Straub, Deshmukh & Elimelech (2016) |
| Gulf seawater salinity | basin mean 40.5 g/kg, coastal values to 45 | Ibrahim et al., Front. Mar. Sci. 7 (2020) 573, https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2020.00573/full |
| Taweelah capacity and energy use | 909,200 m3/day, below 3 kWh/m3 | https://acwapower.com/en/what-we-do/projects/taweelah-ro-desalination-iwp/ and https://www.waterworld.com/drinking-water-treatment/article/14039734/ |
| Fukuoka plant | 10,000 m3/day brine at about 8 %, 9,200 m3/day treated sewage, 110 kW net, 880,000 kWh/year | https://www.f-suiki.or.jp/wp-content/uploads/2023/10/5683f686954e8bd95ac993ef013f1ccb.pdf |
| Abu Dhabi effluent | over 1 million m3/day treated, about 80 % reused, 100 % target | https://www.wam.ae/en/article/blywkjs-abu-dhabi-reuse-100-treated-water-2026-says-taqa |

## Costs and fouling, report sections 8 and 9 (`scripts/s06_economics.py`, `scripts/s07_fouling.py`)

These are first estimates, added after the first draft of the report. The cost figures come from published summaries
and should be checked against the full papers before they are relied on.

| Parameter | Range (base) | Source |
|---|---|---|
| Membrane price | US$5-30/m2 (10) | IRENA, Salinity gradient energy technology brief (2014): EUR 10-30/m2 today, 2-5 needed to compete, https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2014/Salinity_Energy_v4_WEB.pdf ; vendor quotes of US$9.2-12.3/m2 in a 2024 spiral-wound PRO cost study, https://www.sciencedirect.com/science/article/pii/S1944398624202423 |
| Membrane share of PRO capital | 50-80 % (65) | IRENA (2014): membranes are 50-80 % of total capital |
| Membrane life | 3-7 years (5) | 5 years in the 2024 study above; the range is our choice |
| Effluent pretreatment capital | US$150-1000 per m3/day (400) | Order of magnitude only: about US$10 million for a 10,000 m3/day low-pressure membrane plant in 2013, https://dspace.lib.cranfield.ac.uk/server/api/core/bitstreams/f85b9f9c-fc67-4e29-bfb5-ab6a9998ead2/content ; larger plants and shared reuse infrastructure cost less. Weakest cost input. |
| Discount rate, O&M, availability, plant life | 4-10 % (7), 2-4 %/yr (3), 85-95 % (90), 20 years | Our choice |
| Solar comparison | US$0.0135/kWh | Al Dhafra PV bid tariff, Abu Dhabi Power Corporation (2020), https://www.wam.ae/en/article/hszra9cy-abu-dhabi-power-corporation-announces-worlds |
| Published PRO cost check | US$0.44/kWh minimum for RO brine and wastewater | Chung, Banchik, Swaminathan & Lienhard, Desalination 408 (2017) 133-144, https://tethys-engineering.pnnl.gov/publications/present-future-economic-viability-stand-alone-pressure-retarded-osmosis |
| Fouling severity | 0-40 % loss of water permeability, with and without a matching rise in the structural parameter | Our choice, illustrative. No PRO pilot data on Gulf effluent was found; bench forward-osmosis tests on treated wastewater report flux declines of roughly 0-25 %, which the range covers with margin. |

## Values we chose, not measured

The membrane parameters of the lab thin-film composite (A 2.72, B 0.995, S 459) are our own fit to
`straub2014_pro_flux.csv` (see `scripts/s02_validation.py`). The base values and ranges of the plant
inputs are in `scripts/common.py`; the report's parameter table states which are assumptions.
