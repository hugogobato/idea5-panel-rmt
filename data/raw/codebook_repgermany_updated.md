# Codebook: repgermany_updated.dta

## Dataset Information

- **Title**: Updated Replication Data for "Comparative Politics and the Synthetic Control Method"
- **File name**: `repgermany_updated.dta`
- **Format**: Stata (.dta)
- **Dimensions**: 748 observations (17 countries x 44 years), 13 variables
- **Unit of observation**: Country-year
- **Time coverage**: 1960–2003
- **Citation**: Abadie, A., Diamond, A., and Hainmueller, J. (2015). "Comparative Politics and the Synthetic Control Method." *American Journal of Political Science*, 59(2), 495–510.

## Countries

The dataset covers West Germany (the treated unit) and 16 OECD donor countries:

| Index | Country       | Index | Country       |
|-------|---------------|-------|---------------|
| 1     | USA           | 10    | Norway        |
| 2     | UK            | 11    | Switzerland   |
| 3     | Austria       | 12    | Japan         |
| 4     | Belgium       | 13    | Greece        |
| 5     | Denmark       | 14    | Portugal      |
| 6     | France        | 15    | Spain         |
| 7     | West Germany  | 16    | Australia     |
| 8     | Italy         | 17    | New Zealand   |
| 9     | Netherlands   |       |               |

## Variable Descriptions

| Variable | Type | Label | Description | Range | Missing |
|----------|------|-------|-------------|-------|---------|
| `index` | Numeric | Country index | Unique numeric identifier for each country (see table above). | 1–17 | 0 |
| `country` | String | Country name | Name of the country. | — | 0 |
| `year` | Numeric | Year | Calendar year of observation. | 1960–2003 | 0 |
| `gdp` | Numeric | GDP per capita | Per-capita gross domestic product, PPP, current international USD. Source: OECD National Accounts (via OECD Health Database); West Germany data from Statistisches Bundesamt 2005, converted using PPP factors from the OECD Health Database. | 707–37,548 | 0 |
| `infrate` | Numeric | Inflation rate | Annual percentage change in consumer prices (base year 1995). Source: World Bank, World Development Indicators Database 2005. | -0.9–28.8 | 21 |
| `trade` | Numeric | Trade openness | Exports plus imports as percentage of GDP. Source: World Bank, World Development Indicators CD-ROM 2000. | 9.4–149.7 | 102 |
| `schooling` | Numeric | Schooling | Percentage of secondary school attained in total population aged 25+, reported at 5-year intervals. Source: Barro and Lee (2000), CID Working Paper No. 42 — Human Capital Updated Files. | 3.5–69.6 | 597 |
| `invest60` | Numeric | Investment rate (1960s) | Ratio of real domestic investment (private + public) to real GDP, 1960s average. One value per country. Source: Barro and Lee (1994), "Data Set for a Panel of 138 Countries." | 20.5–36.6 | 731 |
| `invest70` | Numeric | Investment rate (1970s) | Ratio of real domestic investment (private + public) to real GDP, 1970s average. One value per country. Source: Barro and Lee (1994), "Data Set for a Panel of 138 Countries." | 21.1–36.8 | 731 |
| `invest80` | Numeric | Investment rate (1980s) | Ratio of real domestic investment (private + public) to real GDP, 1980s average. One value per country. Source: Barro and Lee (1994), "Data Set for a Panel of 138 Countries." | 17.9–31.8 | 731 |
| `industry` | Numeric | Industry share | Industry share of value added (% of GDP). Source: World Bank, World Development Indicators Database 2005. | 21.6–48.0 | 207 |
| `pop` | Numeric | Population | Total population (in thousands). Source: OECD Health Data 2006. | 2,377–290,789 | 0 |
| `deflator_fred` | Numeric | GDP deflator (FRED) | U.S. GDP implicit price deflator (index). Only available for the USA (index = 1). Used for converting current to real (1990) USD. Source: U.S. Bureau of Economic Analysis via FRED (series A191RD3A086NBEA). | 22.2–105.0 | 704 |

## Notes on Missing Values

- `invest60`, `invest70`, `invest80`: These variables contain a single period-average value per country, stored at one year only, with all other country-year observations set to missing.
- `schooling`: Recorded at approximately 5-year intervals only.
- `deflator_fred`: Available only for the USA (44 observations); missing for all other countries.
- `infrate`, `trade`, `industry`: Occasional missing values for some country-years.


