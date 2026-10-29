# wasteWaterVirus

Time-series analysis of virus loads in wastewater and reported cases at two treatment plants in Saxony, Germany (Leipzig and Freiberg).

## Question

Do virus concentrations in sewage track the number of reported cases, and do they lead them in time? The data are weekly virus copy numbers in pooled sewage samples and officially reported cases for SARS-CoV-2, influenza A, astrovirus, rotavirus and norovirus (genotypes I and II), from calendar week 26 of 2021 to week 50 of 2022. The analysis tests whether the two series are correlated, at which lags, and whether adding the sewage signal improves short-term forecasts of cases.

## My role

I designed and ran the time-series analysis for a colleague's wastewater surveillance study: stationarity checks, ARIMA and vector autoregressive (VAR) models, five-week forecasts, and cross-correlation between sewage copy numbers and reported cases. I did not collect or produce the data.

## Status

Not published. The analysis was done for a manuscript draft; the citation will be added here if it is published.

## What the script does

`code.R`, from top to bottom:

1. Reads the ten data files in `data/raw/` and fills missing weeks by linear interpolation (`zoo::na.approx`).
2. For each virus and plant, checks stationarity (ACF, PACF, augmented Dickey-Fuller test, `forecast::ndiffs`) and differences the series where needed.
3. Fits ARIMA models (`forecast::auto.arima`) to copy numbers and to cases separately.
4. Fits VAR models (`vars::VAR`, up to six weekly lags) to copy numbers and cases together and forecasts five weeks ahead.
5. Computes cross-correlations between the differenced copy-number and case series and plots them together, with 95% confidence limits.

## How to run

Last run in October 2026 with R 4.6.0, tidyverse 2.0.0, readxl 1.4.5, forecast 9.0.2, zoo 1.8.15, tseries 0.10.63, lubridate 1.9.5 and vars 1.6.1.

```r
install.packages(c("tidyverse", "readxl", "forecast", "zoo", "tseries", "lubridate", "vars", "bsts"))
```

`bsts` is loaded at the top of the script but not used; install it or comment out that line.

1. Put the ten `.xlsx` files in `data/raw/` (see Data availability). `data/README.txt` describes their columns.
2. Lines 5 and 18 of `code.R` set the working directory and the data path to `~/GitHub/wasteWaterVirus`. Change both to the folder you cloned into.
3. Run `code.R` interactively, top to bottom. It prints the model summaries and draws the diagnostic, forecast and cross-correlation plots; the combined cross-correlation figure is built at the end.

Known issue: the norovirus VAR model with all genotypes (`var6.all`, around line 312) cannot be fitted, because the total norovirus series is the sum of genotypes I and II. `summary(var6.all)` stops with "system is computationally singular". Skip that line and continue.

## Data availability

The data are not included in this repository. They were produced by [colleague's name and institution]; please contact them for access. `data/README.txt` describes the file structure the script expects.

## License

The code is released under the MIT License (see `LICENSE`). The license does not cover the data.
