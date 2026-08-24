# msarima 0.0.0.9006

* New `sacf()` and `spacf()`, autocorrelation and partial autocorrelation plots
  that highlight the lags of each seasonal cycle, with the periods taken from an
  `msts` object or from a fitted model. `ggsacf()` and `ggspacf()` are the
  ggplot2 versions.

* New vignette `msarima`, a walk-through of the package on `MAD_m40`.

* New dataset `MAD_m40`: hourly traffic volume on the M-40 ring road of Madrid,
  February to May 2024, as an `msts` object with seasonal periods 24 and 168.

* Fixed estimation of models with differencing in `msarima()`. The differencing
  operator was folded into the AR polynomial and the conditional sum of squares
  was computed on the undifferenced series, padded with a constant. The first
  `d + sum(D * S)` residuals were therefore start-up transient (they entered the
  likelihood and biased the estimates), and `residuals()` had the same length as
  the input series. The series is now differenced before estimation, so fitting
  `msarima(x, order = c(p, 0, q), seasonal = list(c(0, 1, 1)))` gives exactly the
  same coefficients and residuals as fitting the same ARMA on `diff(x, S)`.
  `residuals()` and `fitted()` are now `NA` for the first `d + sum(D * S)`
  observations (fitted values remain on the scale of the original series).
* `forecast.msarima()` and `simulate.msarima()` now use the `intercept` term of
  models fitted with `include.constant = TRUE`. It was estimated but never read,
  so point forecasts decayed to 0 instead of to the mean of the process, and
  simulated series were centred on 0.
* Fixed seasonal AR/MA polynomial construction in `forecast.msarima()` (seasonal
  orders >= 2 placed the coefficients at lags `S, S+1, ...` instead of
  `S, 2S, ...`) and seasonal differencing for `D >= 2`. The same fix had already
  been applied to `msarima()` and `msarima.sim()`.
* New `intercept` argument in `msarima.sim()`.
* Fixed but to order greater then 1 in `msarima()`
* Fixed bug to extract the hessian in `msarima()`
* Fixed seasonal AR/MA polynomial construction in `msarima()` and
  `msarima.sim()` for seasonal orders >= 2.
* Fixed seasonal AR/MA polynomial construction (and seasonal
  differencing for D >= 2) in `msarima()` and `msarima.sim()`.
* Initial CRAN submission.
