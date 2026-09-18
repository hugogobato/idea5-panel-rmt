# Dataset Scout Report — Fragile / Sub-frontier SCM–SDID Panels

**Context:** BBP-type spectral recoverability frontier (m = s/√c = 1) for synthetic-control / matrix-completion estimators on low-rank factor panels Y = L + E. Pipeline: outcome-only, unit-centered levels, n_d ≥ ~6 donors, T0 ≥ ~12 pre-periods. Canonical macro showcases (CA Prop 99, German reunification, Basque) saturate the frontier because shared stochastic trends make macro levels supercritical BY CONSTRUCTION — this list targets panels **not** dominated by a common trend, with documented or contested fit.

**Ranked candidates** (rank 1 = most promising fragile case). "†" = claim not fully re-verifiable from the fetched source this session (see Honesty flags).

| Rank | Paper (year, journal) | Treated unit / Treatment | Outcome | n_d × T0 (c = n_d/T0) | Archive DOI/URL + license | Why plausibly sub-frontier / misaligned (fit quote) | Access notes |
|---|---|---|---|---|---|---|---|
| 1 | **Powell (2026), "Imperfect Synthetic Controls"** (PMC open access) [1] | Wisconsin — repeal of 48-h handgun purchase waiting period (Jun 2015) | Handgun-suicide deaths per 100,000 (monthly) | ~49 other states × ~161 monthly pre-periods (~13 yr) → **c ≈ 0.30** ✓ | [PMC12843250](https://pmc.ncbi.nlm.nih.gov/articles/PMC12843250/) — open access (CC-BY-NC) | Documented poor pre-fit + convex-hull violation; micro outcome, no common-trend saturation. *"While Wisconsin's synthetic control is a poor fit…"* [1] | No packaged .csv verified; reconstruct from CDC WONDER (public). Single treated unit, long clean pre-window. |
| 2 | **Donohue, Aneja & Weber (2019), "Right-to-Carry Laws and Violent Crime"** (J. Empirical Legal Studies) [2] | 33 states adopting Right-to-Carry (staggered, ~1977–2010s) | Violent crime & murder rates | ~30 never-adopting states × up to ~28 pre-yr (early adopters) → **c ≈ 1** ✓ | [NBER w23510](https://www.nber.org/system/files/working_papers/w23510/revisions/w23510.rev0.pdf); packaged archive unverified | Contested fit — Moody & Marvell (2019) [3]: *"the synthetic control model fails to control for any of the major factors that cause crime rates to vary"*; no net effect across all 33 states. | Slice early-adopter states into single-treated runs; donors = never-adopting states. Crime = micro outcome. Also in staggered list. |
| 3 | **Moser & Voena (2012), "Compulsory Licensing: Evidence from the Trading with the Enemy Act"** (AER) [4] | Germany — 1919 compulsory licensing of US/Allied patents | Domestic invention (patent counts) | ~20–30 countries × ~15 pre-yr (unverified) | [openICPSR 112497](https://www.openicpsr.org/openicpsr/project/112497/version/V1/view) — `chem_patents_maindataset.dta` confirmed; public | Patent counts are noisy micro data, not a smooth common trend; plausibly sub-frontier. Fit quote: not retrieved † | Archive confirmed with .dta; reshape country×year. Treatment timing (1919) unambiguous. |
| 4 | **Jardim et al. (2022), "Minimum Wage Increases and Low-Wage Employment: Evidence from Seattle"** (AEJ: Policy) [5] | Seattle — 2015/2016 wage hikes | Hours worked & wages (micro) | ~100+ comparison PUMAs × ~20 quarterly pre → **c ≈ 5** (unverified) | AEJ: Policy replication (AEA); file list unverified † | Micro labor outcome; heavily contested by replications. Fit quote: not retrieved † | Quarterly panel likely clears T0 ≥ 12 — verify periodicity before running. |
| 5 | **Skorobogatov (2021), "The effect of alcohol sales restrictions on alcohol poisoning mortality: Evidence from Russia"** (Health Economics) [6] | 9 Russian regions (restrictions by 2005/06) | log alcohol-poisoning mortality | 47 donors (13 pos-weight) × **9 pre-yr** → c = 5.2 (**T0 below floor**) | No archive — Rosstat (public) + Consultant Plus (commercial) | Explicitly documented poor pre-fit. *"This is mostly explained by the poor fit of alcohol poisoning mortality by the SC prior to the treatment period for the respective regions from the control group"*; also *"several regions with a relatively poor fit… explained by their unusual levels of alcohol poisoning mortality which cannot be reproduced."* [6] | Reconstruct from Rosstat; treatment dates in paper. Needs longer monthly variant or flag. |
| 6 | **Gobillon & Magnac (2016), "Regional Policy Evaluation: Interactive Fixed Effects and Synthetic Controls"** (REStat) [7] | 13 Paris-region municipalities (ZFU enterprise zones, created Jan 1997) | log unemployment-exit ratio (semester) | 135 control municipalities × **8 pre-semesters** → c = 16.9 (**T0 below floor**) | [TSE WP PDF](https://publications.ut-capitole.fr/15682/1/wp_tse_419.pdf); packaged archive unverified | Methods paper WITH real application; authors distrust SCM here. *"synthetic controls estimates are more sensitive to the specification than factor models and difference-in-differences estimates… we have reasons to believe that interactive effect estimates are more credible."* [7] | Micro outcome (employment exits); French admin data reconstruction needed. |
| 7 | **Peri & Yasenov (2019), "The Labor Market Effects of a Refugee Wave: SCM Meets the Mariel Boatlift"** (J. Human Resources) [8] | Miami — 1980 Mariel Boatlift | log weekly/hourly wages, unemployment | 43 cities (May-ORG) / 31 (March-CPS) × **6 pre-yr** → c = 7.2 (**T0 below floor**) | [JHR article](https://jhr.uwpress.org/content/54/2/267) + [author appendix PDF](https://giovanniperi.ucdavis.edu/uploads/5/6/8/2/56826033/mariel_jun2_app.pdf) — public | Most famous contested-fit case (Borjas 2017 vs PY); sensitivity to matching window is the crux. *"Most results are robust to variations in the selection of the Synthetic Control Group, as long as one matches the whole 1973…"* [8] | Monthly CPS variants can extend T0; annual main spec is below floor. |
| 8 | **Dube & Zipperer (2015), "Pooling Multiple Case Studies Using Synthetic Controls"** (IZA DP 8944) [9] | 29 state minimum-wage cases (1979–2013) | Teen wage & employment | 29 cases × case-specific T0 | [arindube.com replication](https://arindube.com/working-papers/) — aggregated panel, public | Pooled SCM with documented cross-case heterogeneity; good for donor-thinning / many-case robustness. Fit quote: not retrieved † | Archive confirmed present (aggregated panel). |
| 9 | **Fletcher, Frisvold & Tefft (2015), "Non-linear effects of soda taxes on consumption and weight outcomes"** (Health Economics) [10] | Soda-tax states | Soda consumption & weight/BMI | ~40 states × **4 pre-yr** → c ≈ 10 (**T0 below floor**) | No packaged archive verified; public survey data | Documented poor fit for a treated state. *"is significant only because the overall synthetic control 'fit' is poor for Ohio."* [10] | T0 = 4 far below floor — only usable with a longer monthly/quarterly variant. |
| 10 | **Courtemanche & Zapata (2014), "Does Universal Coverage Improve Health? The Massachusetts Experience"** (JPAM) [11] | Massachusetts — 2006 health reform | Mortality | ~40 states × ~5 pre-yr → c ≈ 8 (**T0 below floor**) | [NBER w17893](https://www.nber.org/system/files/working_papers/w17893/w17893.pdf) — public WP; packaged archive unverified | Contested fit — Powell (2026) [1] critique of this application: *"donor states had poor pretreatment match on predictor variables."* | Mortality micro outcome; short T0. |
| 11 | **Cavallo, Galiani, Noy & Pantano (2013), "Catastrophic Natural Disasters and Economic Growth"** (REStat) [12] | Disaster countries (e.g., Chile 1960) | GDP **growth** (not levels) | ~20 countries × ~20 yr (unverified) | [IDB publication](https://publications.iadb.org/en/catastrophic-natural-disasters-and-economic-growth); PWT data (public) | Headline result is "no effect" (weak signal); growth rates avoid the world-cycle level saturation. Dimensions & fit not verified † | Macro, but growth not levels — lower priority. |

## Top-3 fit-for-pipeline

1. **Wisconsin handgun suicide (Powell 2026)** — cleanest: single treated unit, ~161 monthly pre-periods (comfortably clears the T0 ≥ 12 floor), ~49 donors, micro outcome, poor fit documented in-paper. Reshaping: Wisconsin vs other-states monthly matrix from CDC WONDER; treatment date (Jun 2015) unambiguous.
2. **RTC laws (DAW 2019)** — staggered set where early-adopter states give long pre-windows (T0 ≈ 20–28 yr) as separate single-treated runs with never-adopting states as donors — the ideal c ≈ 1 regime. Contested fit (Moody–Marvell [3]) is exactly the misaligned signal wanted; main cost is assembling the UCR state panel (no verified packaged archive).
3. **Moser–Voena (2012)** — best *archived* micro-outcome panel (verified `chem_patents_maindataset.dta` on openICPSR [4]); single treated unit (Germany, 1919); noisy patent counts rather than smooth trend. Fit unverified † — run a quick BBP check first; if supercritical, demote to negative control.

## Negative controls (expected clearly supercritical)

1. **California Prop 99 smoking** — Abadie, Diamond & Hainmueller (2010, JASA) [13]; good fit; already supercritical in project runs.
2. **German reunification** — ADH (2015, AJPS) [14], Harvard Dataverse 24714; good fit; already supercritical.
3. **Economic liberalization episodes** — Billmeier & Nannicini (2013, REStat) [15]; national GDP-per-capita *levels* with a dominant world-cycle factor — the exact "supercritical by construction" failure mode already saturated.

## Staggered-adoption candidates (fragile profile)

- **Right-to-Carry laws** (DAW 2019) [2] — 33 states, contested fit (top staggered pick).
- **Medicaid expansion & mortality** (Courtemanche et al. 2023, NBER w30818) [16] — staggered states, mortality, generalized SCM; URL unverified †.
- **Minimum wage** (Dube–Zipperer 2015) [9] — 29 staggered cases, teen wage/employment.
- **Marijuana legalization & crime** (Harper 2023, J. Drug Issues) [17] — staggered states, crime rates.
- **Naloxone co-prescribing laws** (Duska et al. 2022) [18] — 5 states, opioid-overdose mortality; full citation unverified †.

## Honesty flags (unverified this session)

- **Fit quotes** for Moser–Voena [4], Seattle/Jardim [5], Dube–Zipperer [9], Harper [17]: none retrieved — marked †.
- **Wisconsin specifics** (161-month figure, convex-hull sentence, 49-donor count): the PMC fetch was truncated, so only the snippet-level quote *"While Wisconsin's synthetic control is a poor fit…"* is verbatim-confirmed; re-verify the rest against the PDF before citing.
- **Packaged rectangular archives confirmed**: Moser–Voena (openICPSR `.dta` [4]), Dube–Zipperer (arindube.com [9]), Peri–Yasenov (JHR/author page [8]). **Not yet confirmed** (crawler blocked): openICPSR 117261 (Cherry Picking [19], project page verified via search), openICPSR 146381 (SDID [20] — project ID unverified), AEA article pages, RTC, Seattle, MA-health, Skorobogatov, Gobillon–Magnac, Fletcher, Cavallo.
- **Pipeline floor (n_d ≥ 6, T0 ≥ 12)**: only Wisconsin, RTC (early adopters), Seattle (quarterly), and likely Moser–Voena clear it; all others flagged per-row.

## References

[1] Powell, D. (2026). "Imperfect Synthetic Controls." NIH/PMC open access. https://pmc.ncbi.nlm.nih.gov/articles/PMC12843250/

[2] Donohue, J. J., III, Aneja, A., & Weber, K. D. (2019). "Right-to-Carry Laws and Violent Crime: A Comprehensive Assessment Using Panel Data and a State-Level Synthetic Controls Analysis." *Journal of Empirical Legal Studies*. NBER WP version: https://www.nber.org/system/files/working_papers/w23510/revisions/w23510.rev0.pdf

[3] Moody, C. E., & Marvell, T. B. (2019). "Do Right-to-Carry Laws Increase Violent Crime? A Comment on Donohue, Aneja, and Weber." *Econ Journal Watch*. https://econjwatch.org/articles/do-right-to-carry-laws-increase-violent-crime-a-comment-on-donohue-aneja-and-weber

[4] Moser, P., & Voena, A. (2012). "Compulsory Licensing: Evidence from the Trading with the Enemy Act." *American Economic Review*, 102(1), 396–427. AEA: https://www.aeaweb.org/articles?id=10.1257/aer.102.1.396 — Replication (openICPSR 112497, file `chem_patents_maindataset.dta`): https://www.openicpsr.org/openicpsr/project/112497/version/V1/view

[5] Jardim, E., Long, M. C., Plotnick, R., van Inwegen, E., Vigdor, J., & Wething, H. (2022). "Minimum Wage Increases and Low-Wage Employment: Evidence from Seattle." *AEJ: Economic Policy*. (Article/replication URL not verified this session †)

[6] Skorobogatov, A. S. (2021). "The effect of alcohol sales restrictions on alcohol poisoning mortality: Evidence from Russia." *Health Economics*. Wiley: https://onlinelibrary.wiley.com/doi/abs/10.1002/hec.4251 — Author PDF: https://game.hse.ru/data/2020/03/10/1562997014/skorobogatov_paper.pdf

[7] Gobillon, L., & Magnac, T. (2016). "Regional Policy Evaluation: Interactive Fixed Effects and Synthetic Controls." *Review of Economics and Statistics*, 98(3), 535–551. MIT Press: https://direct.mit.edu/rest/article/98/3/535/58348 — WP PDF: https://publications.ut-capitole.fr/15682/1/wp_tse_419.pdf

[8] Peri, G., & Yasenov, V. (2019). "The Labor Market Effects of a Refugee Wave: Synthetic Control Method Meets the Mariel Boatlift." *Journal of Human Resources*, 54(2), 267. https://jhr.uwpress.org/content/54/2/267 — Appendix PDF: https://giovanniperi.ucdavis.edu/uploads/5/6/8/2/56826033/mariel_jun2_app.pdf

[9] Dube, A., & Zipperer, B. (2015). "Pooling Multiple Case Studies Using Synthetic Controls: An Application to Minimum Wage Policies." IZA DP 8944. https://ideas.repec.org/p/iza/izadps/dp8944.html — Replication data & programs: https://arindube.com/working-papers/

[10] Fletcher, J. M., Frisvold, D. E., & Tefft, N. (2015). "Non-linear effects of soda taxes on consumption and weight outcomes." *Health Economics*. Wiley: https://onlinelibrary.wiley.com/doi/abs/10.1002/hec.3045 — NIHMS PDF: https://pmc.ncbi.nlm.nih.gov/articles/PMC6047515/pdf/nihms979522.pdf

[11] Courtemanche, C. J., & Zapata, D. (2014). "Does Universal Coverage Improve Health? The Massachusetts Experience." *Journal of Policy Analysis and Management*, 33(1), 160–185. NBER WP: https://www.nber.org/system/files/working_papers/w17893/w17893.pdf

[12] Cavallo, E., Galiani, S., Noy, I., & Pantano, J. (2013). "Catastrophic Natural Disasters and Economic Growth." *Review of Economics and Statistics*, 95(5), 1549–1561. IDB mirror: https://publications.iadb.org/en/catastrophic-natural-disasters-and-economic-growth

[13] Abadie, A., Diamond, A., & Hainmueller, J. (2010). "Synthetic Control Methods for Comparative Case Studies: Estimating the Effect of California's Tobacco Control Program." *Journal of the American Statistical Association*, 105(490), 493–505. NBER WP: https://www.nber.org/system/files/working_papers/w12831/revisions/w12831.rev0.pdf

[14] Abadie, A., Diamond, A., & Hainmueller, J. (2015). "Comparative Politics and the Synthetic Control Method." *American Journal of Political Science*, 59(2), 495–510. Replication (Harvard Dataverse, doi:10.7910/DVN/24714): https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/24714

[15] Billmeier, A., & Nannicini, T. (2013). "Assessing Economic Liberalization Episodes: A Synthetic Control Approach." *Review of Economics and Statistics*, 95(3), 983–1001. (URL not verified this session †)

[16] Courtemanche, C., et al. (2023). ACA Medicaid Expansion and mortality. NBER WP 30818. (URL not verified this session †)

[17] Harper, A. J. (2023). "Estimating the Effect of Legalizing Marijuana on Crime Rates." *Journal of Drug Issues*. https://journals.sagepub.com/doi/abs/10.1177/00220426221134107

[18] Duska, M., et al. (2022). "Naloxone co-prescribing laws and opioid overdose mortality." 5-state SCM/SDID application. (Full citation and URL not verified this session †)

[19] Ferman, B., Pinto, C., & Possebom, V. (2020). "Cherry Picking with Synthetic Controls." *Journal of Policy Analysis and Management*, 39(2), 510–532. https://onlinelibrary.wiley.com/doi/abs/10.1002/pam.22206 — Replication (openICPSR 117261): https://www.openicpsr.org/openicpsr/project/117261/version/V2/view

[20] Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021). "Synthetic Difference-in-Differences." *American Economic Review*, 111(12), 4088–4118. https://www.aeaweb.org/articles?id=10.1257/aer.20190159 — Replication (openICPSR; project ID unverified †): https://www.openicpsr.org/openicpsr/project/146381/version/V1/view
