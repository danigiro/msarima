# msarima: Multi-Seasonal ARIMA Model

Fit, forecast and simulate autoregressive integrated moving average
(ARIMA) models with more than one seasonal component, following Lisi and
Grigoletto (2025)
[doi:10.1002/for.3283](https://doi.org/10.1002/for.3283) . Series
observed at high frequency often display several seasonal cycles at
once, such as the daily and the weekly pattern of hourly data, which a
single seasonal ARIMA operator cannot represent. Estimation is by
conditional maximum likelihood through numerical optimisation, with
support for fixed coefficients and a constant term. Seasonal
correlograms, forecasts with prediction intervals, and simulation from a
model specification or from a fitted model are also provided.

## See also

Useful links:

- <https://github.com/danigiro/msarima>

- Report bugs at <https://github.com/danigiro/msarima/issues>

## Author

**Maintainer**: Daniele Girolimetto <daniele.girolimetto@unipd.it>
([ORCID](https://orcid.org/0000-0001-9387-1232)) \[copyright holder\]

Authors:

- Francesco Lisi <francesco.lisi@unipd.it>
  ([ORCID](https://orcid.org/0000-0002-8749-8889))

- Matteo Grigoletto
