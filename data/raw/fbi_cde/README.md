# FBI Crime Data Explorer monthly extract

`fbi_cde_monthly_1999_2017.csv` was acquired on 2026-08-29 from the official FBI Crime Data Explorer web endpoint by running:

```bash
python3 scripts/fetch_fbi_cde.py
```

The query requests all 50 states, `violent-crime`, `01-1999` through `12-2017`, and `type=estimates`. It returns 11,400 state-month rows. The API's nested `actuals` field is retained as `actual`; under `type=estimates`, this is the CDE estimated monthly count, not a claim that every agency reported. `coverage_pct` and `participated_population` are retained so policy panels can be screened without outcome leakage.

The acquisition URL template is:

```text
https://cde.ucr.cjis.gov/LATEST/summarized/state/{STATE}/violent-crime?from=01-1999&to=12-2017&type=estimates
```

No API key was required on the acquisition date. The sidecar metadata file records every resolved request URL and the returned CDE refresh metadata. The loader `applications.loaders.load_fbi_crime_proxy` requires a complete rectangular state-month grid and applies any minimum-coverage rule using pre-treatment months only.

SHA-256 checksums are registered in `data/raw/SHA256SUMS`. FBI UCR participation is voluntary, so coverage filtering and composition sensitivity are required in every applied use. A rectangular panel alone is not a data-quality pass.
