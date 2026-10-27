# Simulate from a Multi-Seasonal ARIMA Model (Direct Interface)

Simulates a time series from a multi-seasonal ARIMA model with specified
parameters. This function provides a direct interface for simulation
without requiring a fitted model object.

## Usage

``` r
msarima.sim(
  nsim,
  order = c(0, 0, 0),
  seasonal = list(),
  seasonal.period = NULL,
  ar = NULL,
  ma = NULL,
  sigma2 = 1,
  burn = NULL,
  intercept = 0
)
```

## Arguments

- nsim:

  Length of the simulated series (required).

- order:

  A numeric vector of length 3 specifying the non-seasonal part of the
  ARIMA model: the AR order (p), the degree of differencing (d), and the
  MA order (q). Default is `c(0, 0, 0)`.

- seasonal:

  A list of numeric vectors, each of length 3, specifying the seasonal
  components. Each vector gives the seasonal AR order (P), the seasonal
  differencing order (D), and the seasonal MA order (Q) of one
  component. Use [`list()`](https://rdrr.io/r/base/list.html) for a
  model without seasonal components. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- seasonal.period:

  A numeric vector of seasonal periods, one for each element of
  `seasonal`. Periods can be given in any order but must be distinct.
  Default is `NULL`.

- ar:

  Numeric vector of AR coefficients, non-seasonal first and then one
  block per seasonal component. Its length must match the total AR order
  implied by `order` and `seasonal`. Default is `NULL`.

- ma:

  Numeric vector of MA coefficients, ordered as `ar`. Its length must
  match the total MA order implied by `order` and `seasonal`. Default is
  `NULL`.

- sigma2:

  Innovation variance. Default is 1.

- burn:

  Number of burn-in observations to discard. If `NULL`, it is set to
  `max(100, 2 * kk)`, where `kk` is the largest lag in the model.
  Default is `NULL`.

- intercept:

  Constant term of the model equation, that is \\\mu\\ in \\y_t = \mu +
  \sum \phi_i y\_{t-i} + \sum \theta_j e\_{t-j} + e_t\\. This is the
  intercept, not the mean of the process: for a stationary model the
  mean is `intercept / prod(1 - ar)`. Default is 0.

## Value

A time series object of length `nsim`: `ts` with one seasonal period,
`msts` with several.

## See also

[`msarima`](https://danigiro.github.io/msarima/reference/msarima.md),
[`simulate.msarima`](https://danigiro.github.io/msarima/reference/simulate.msarima.md)

## Examples

``` r
# Simulate data
set.seed(123)
y <- msarima.sim(nsim = 200,
                 order = c(1, 0, 0),
                 seasonal = list(c(1, 0, 0), c(0, 0, 1)),
                 seasonal.period = c(4, 11),
                 ar = c(0.5, 0.5), ma = 0.3, sigma2 = 1)
```
