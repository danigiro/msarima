# Changelog

## msarima 0.0.0.9006

- New [`sacf()`](https://danigiro.github.io/msarima/reference/sacf.md)
  and [`spacf()`](https://danigiro.github.io/msarima/reference/sacf.md),
  autocorrelation and partial autocorrelation plots that highlight the
  lags of each seasonal cycle, with the periods taken from an `msts`
  object or from a fitted model.
  [`ggsacf()`](https://danigiro.github.io/msarima/reference/sacf.md) and
  [`ggspacf()`](https://danigiro.github.io/msarima/reference/sacf.md)
  are the ggplot2 versions.

- New vignette `msarima`, a walk-through of the package on `MAD_m40`.

- New dataset `MAD_m40`: hourly traffic volume on the M-40 ring road of
  Madrid, February to May 2024, as an `msts` object with seasonal
  periods 24 and 168.

- Fixed estimation of models with differencing in
  [`msarima()`](https://danigiro.github.io/msarima/reference/msarima.md).
  The differencing operator was folded into the AR polynomial and the
  conditional sum of squares was computed on the undifferenced series,
  padded with a constant. The first `d + sum(D * S)` residuals were
  therefore start-up transient (they entered the likelihood and biased
  the estimates), and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) had the same
  length as the input series. The series is now differenced before
  estimation, so fitting
  `msarima(x, order = c(p, 0, q), seasonal = list(c(0, 1, 1)))` gives
  exactly the same coefficients and residuals as fitting the same ARMA
  on `diff(x, S)`.
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) are now `NA`
  for the first `d + sum(D * S)` observations (fitted values remain on
  the scale of the original series).

- [`forecast.msarima()`](https://danigiro.github.io/msarima/reference/forecast.msarima.md)
  and
  [`simulate.msarima()`](https://danigiro.github.io/msarima/reference/simulate.msarima.md)
  now use the `intercept` term of models fitted with
  `include.constant = TRUE`. It was estimated but never read, so point
  forecasts decayed to 0 instead of to the mean of the process, and
  simulated series were centred on 0.

- Fixed seasonal AR/MA polynomial construction in
  [`forecast.msarima()`](https://danigiro.github.io/msarima/reference/forecast.msarima.md)
  (seasonal orders \>= 2 placed the coefficients at lags `S, S+1, ...`
  instead of `S, 2S, ...`) and seasonal differencing for `D >= 2`. The
  same fix had already been applied to
  [`msarima()`](https://danigiro.github.io/msarima/reference/msarima.md)
  and
  [`msarima.sim()`](https://danigiro.github.io/msarima/reference/msarima.sim.md).

- New `intercept` argument in
  [`msarima.sim()`](https://danigiro.github.io/msarima/reference/msarima.sim.md).

- Fixed but to order greater then 1 in
  [`msarima()`](https://danigiro.github.io/msarima/reference/msarima.md)

- Fixed bug to extract the hessian in
  [`msarima()`](https://danigiro.github.io/msarima/reference/msarima.md)

- Fixed seasonal AR/MA polynomial construction in
  [`msarima()`](https://danigiro.github.io/msarima/reference/msarima.md)
  and
  [`msarima.sim()`](https://danigiro.github.io/msarima/reference/msarima.sim.md)
  for seasonal orders \>= 2.

- Fixed seasonal AR/MA polynomial construction (and seasonal
  differencing for D \>= 2) in
  [`msarima()`](https://danigiro.github.io/msarima/reference/msarima.md)
  and
  [`msarima.sim()`](https://danigiro.github.io/msarima/reference/msarima.sim.md).

- Initial CRAN submission.
