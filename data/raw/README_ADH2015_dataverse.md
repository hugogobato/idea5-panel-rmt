# Replication Archive

## Erratum: Comparative Politics and the Synthetic Control Method

**Alberto Abadie**, Department of Economics, MIT
**Alexis Diamond**, College of Computational Sciences, Minerva University
**Jens Hainmueller**, Department of Political Science, Stanford University

February 2026

This archive contains updated replication code and data for:

> Abadie, A., Diamond, A., and Hainmueller, J. (2015). "Comparative Politics and the Synthetic Control Method." *American Journal of Political Science*, 59(2), 495–510.

## Summary of Corrections

1. **Outcome label**: The original article reported GDP in PPP-adjusted 2002 USD, whereas the analysis used GDP in PPP-adjusted current USD. The mislabeling reflects language carried over from an earlier draft that used a different dataset.

2. **OECD predictor averages**: A coding error affected the population-weighted averages in the last column of Tables 2 and 4. The correction affects only these descriptive statistics and no other results.

## Instructions for Replication

1. Ensure all required software and packages are installed (see below).
2. Set the R working directory to the root of this replication archive (the folder containing this README and `repgermany_updated.dta`).
3. Run the three R scripts in any order. Each script automatically creates its output directory if it does not already exist.

| Script | Description | Output Directory | Approximate Runtime |
|--------|-------------|-----------------|---------------------|
| `rep_updated.R` | Main replication (Appendix Section I: Tables 1–4, Figures 1–7) with corrected labeling and OECD averages. | `results/` | ~3 minutes |
| `rep_updated_real_PPP_ex_post.R` | Appendix Section II: Ex post deflated per-capita GDP (Tables A1–A4, Figures A1–A7). Deflation applied after constructing the synthetic control. | `results_real_ex_post/` | ~3 minutes |
| `rep_updated_real_PPP_ex_ante.R` | Appendix Section III: Ex ante deflated per-capita GDP (Tables B1–B4, Figures B1–B7). GDP deflated before constructing the synthetic control. | `results_real_ex_ante/` | ~3 minutes |

Total runtime for all three scripts is approximately 10 minutes.

## Computing Environment

The code was tested in the following environment:

- **Operating System**: macOS (Darwin), Apple M2 Ultra chip
- **Processor**: Apple M2 Ultra, 24 cores (16 Performance + 8 Efficiency)
- **Memory (RAM)**: 192 GB
- **Software**: R 4.5.2 (2025-10-31), platform aarch64-apple-darwin20

## Software Requirements

- R (version 4.5.x or later recommended)
- R packages (with versions used for testing):
  - `foreign` (0.8-91)
  - `Synth` (1.1-9)
  - `xtable` (1.8-8)
  - `dplyr` (1.2.0)
  - `gtools` (3.9.5)
  - `kernlab` (0.9-33)

## Input Data

| File | Description |
|------|-------------|
| `repgermany_updated.dta` | Panel dataset with per-capita GDP (PPP, current USD) and covariates for West Germany and 16 OECD donor countries, 1960–2003. Includes population and U.S. GDP deflator for real GDP calculations. See `codebook_repgermany_updated.md` for variable definitions. |

## Data Sources

The analytic dataset `repgermany_updated.dta` is based on the original replication dataset from Abadie, Diamond, and Hainmueller (2015), supplemented with additional variables (population, GDP deflator) for the corrections in this erratum.

1. **Original replication data**:

   Abadie, A., Diamond, A., and Hainmueller, J. (2015). "Comparative Politics and the Synthetic Control Method." *American Journal of Political Science*, 59(2), 495–510. Replication data available at: https://doi.org/10.7910/DVN/24714, Harvard Dataverse.

   The original dataset draws on the following underlying sources:

   - **GDP per capita** (PPP, current USD): OECD National Accounts (retrieved via the OECD Health Database). Data for West Germany obtained from Statistisches Bundesamt 2005 (Arbeitskreis "Volkswirtschaftliche Gesamtrechnungen der Lander") and converted using PPP monetary conversion factors (retrieved from the OECD Health Database).
   - **Trade openness** (exports plus imports as % of GDP): World Bank, *World Development Indicators* CD-ROM 2000.
   - **Inflation rate** (annual % change in consumer prices, base year 1995): World Bank, *World Development Indicators* Database 2005.
   - **Industry share** (industry share of value added, % of GDP): World Bank, *World Development Indicators* Database 2005.
   - **Schooling** (% of secondary school attained in total population aged 25+, reported in 5-year increments): Barro, R. J., and Lee, J.-W. (2000). "International Data on Educational Attainment: Updates and Implications." CID Working Paper No. 42 — Human Capital Updated Files.
   - **Investment rate** (ratio of real domestic investment to real GDP, reported in 5-year averages): Barro, R. J., and Lee, J.-W. (1994). "Data Set for a Panel of 138 Countries." Available at http://www.nber.org/pub/barro.lee/.

2. **Population data**:

   OECD Health Data 2006. Organisation for Economic Co-operation and Development. Total population for OECD countries, 1960–2003.

3. **U.S. GDP deflator** (used for real GDP calculations in Appendix Sections II and III):

   U.S. Bureau of Economic Analysis. "Gross Domestic Product: Implicit Price Deflator" [A191RD3A086NBEA]. Retrieved from FRED, Federal Reserve Bank of St. Louis. https://fred.stlouisfed.org/series/A191RD3A086NBEA.
