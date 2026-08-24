# Simulate from a Multi-Seasonal ARIMA Model

Returns a simulated time series from a fitted ARIMA model with one or
more seasonal components.

## Usage

``` r
# S3 method for class 'msarima'
simulate(object = NULL, nsim = 100, seed = NULL, burn = 100, ...)
```

## Arguments

- object:

  An object of class `"msarima"`.

- nsim:

  Length of the simulated series. Default is 100.

- seed:

  Either `NULL` or an integer passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html) before simulating, to
  make the result reproducible. Default is `NULL`.

- burn:

  Number of burn-in observations to discard. Default is 100.

- ...:

  Additional arguments (currently not used).

## Value

A time series object (`ts` or `msts`) of length `nsim`, with the
seasonal periods of the fitted model.

## Details

The model specification and the estimated parameters, including the
constant term when the model was fitted with `include.constant = TRUE`,
are extracted from the fitted object and passed to
[`msarima.sim`](https://danigiro.github.io/msarima/reference/msarima.sim.md).

## See also

[`msarima`](https://danigiro.github.io/msarima/reference/msarima.md),
[`msarima.sim`](https://danigiro.github.io/msarima/reference/msarima.sim.md)

## Examples

``` r
# Simulate data
set.seed(123)
y <- msarima.sim(nsim = 200,
                 order = c(1, 0, 0),
                 seasonal = list(c(1, 0, 0), c(0, 0, 1)),
                 seasonal.period = c(4, 11),
                 ar = c(0.5, 0.5), ma = 0.3, sigma2 = 1)

# Fit model
fit <- msarima(y, order = c(1, 0, 0),
               seasonal = list(c(1, 0, 0), c(0, 0, 1)),
               seasonal.period = c(4, 11))

# Simulate from fitted model
y_sim <- simulate(fit, nsim = 100)
```
