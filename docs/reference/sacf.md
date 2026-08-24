# Autocorrelation Plots for Multiple Seasonal Cycles

Sample autocorrelation (`sacf`) and partial autocorrelation (`spacf`)
functions that highlight the lags belonging to one or more seasonal
cycles, so that the seasonal structure of a multi-seasonal series is
visible at a glance. `ggsacf()` and `ggspacf()` are the ggplot2
equivalents.

## Usage

``` r
sacf(
  x,
  season = NULL,
  lag.max = NULL,
  ci = 0.95,
  col = NULL,
  main = NULL,
  ylab = "ACF",
  plot = TRUE,
  ...
)

spacf(x, ..., ylab = "PACF")

ggsacf(
  x,
  season = NULL,
  lag.max = NULL,
  ci = 0.95,
  main = NULL,
  ylab = "ACF",
  ...
)

ggspacf(x, ..., ylab = "PACF")
```

## Arguments

- x:

  A numeric vector, a `ts`/`msts` object, or a fitted
  [msarima](https://danigiro.github.io/msarima/reference/msarima.md)
  model (in which case the model residuals are used).

- season:

  Numeric vector of seasonal periods to highlight, of any length. If
  `NULL` (default) the periods are taken from the `msts`/`ts` attributes
  of `x`, or from the fitted model.

- lag.max:

  Maximum lag. Defaults to `3 * max(season)`, or to the
  [stats::acf](https://rdrr.io/r/stats/acf.html) default when no
  seasonal period is available.

- ci:

  Coverage of the confidence bands, as a probability. Default `0.95`.

- col:

  Colours used for the seasonal lags, one per element of `season`
  (recycled). Defaults to a colour-blind safe palette (Okabe and Ito,
  2008), the same one used by `ggsacf()` and `ggspacf()`.

- main, ylab:

  Plot title and y-axis label.

- plot:

  If `TRUE` (default) the correlogram is drawn.

- ...:

  Further arguments passed to
  [stats::acf](https://rdrr.io/r/stats/acf.html), and from `spacf()` and
  `ggspacf()` to `sacf()` and `ggsacf()`.

## Value

`sacf()` and `spacf()` invisibly return the
[stats::acf](https://rdrr.io/r/stats/acf.html) object, with the seasonal
periods in the `"season"` attribute and, for each lag, the index of the
period it belongs to in `"season.group"`. `ggsacf()` and `ggspacf()`
return a `ggplot` object.

## Details

Lags that are multiples of more than one period (lag 168 with
`season = c(24, 168)`, say) are attributed to the longest matching
period.

The series is treated as a plain numeric vector, so lags are always
reported in number of observations and not in units of time. Leading and
trailing `NA`s, such as those produced by differencing in
[msarima](https://danigiro.github.io/msarima/reference/msarima.md), are
dropped
([stats::na.contiguous](https://rdrr.io/r/stats/na.contiguous.html)).

## See also

[stats::acf](https://rdrr.io/r/stats/acf.html),
[msarima](https://danigiro.github.io/msarima/reference/msarima.md)

## Examples

``` r
fit <- msarima(MAD_m40, order = c(1, 0, 2),
               seasonal = list(c(0, 1, 1), c(0, 1, 1)))

sacf(MAD_m40, lag.max = 200)   # periods 24 and 168 taken from the msts object

spacf(fit)                     # residual pacf, periods taken from the model

```
