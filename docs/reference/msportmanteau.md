# Multi-Seasonal Portmanteau Tests

Multi-seasonal generalisation of the Ljung-Box and Monti portmanteau
tests: residual autocorrelation is checked at low-order and seasonal
lags at once, for an arbitrary number of seasonal periods.

## Usage

``` r
msportmanteau(
  x,
  type = c("Ljung-Box", "Monti"),
  season = NULL,
  nlag = 5,
  lags_short = 0,
  npar = NULL
)
```

## Arguments

- x:

  A numeric vector, a `ts`/`msts` object, or a fitted
  [msarima](https://danigiro.github.io/msarima/reference/msarima.md)
  model (in which case the model residuals are used).

- type:

  Either `"Ljung-Box"` (default), which builds the statistic from the
  autocorrelations, or `"Monti"`, which uses the partial
  autocorrelations. Partial matching is allowed.

- season:

  Numeric vector of seasonal periods to test, of any length. If `NULL`
  (default) the periods are taken from the `msts`/`ts` attributes of
  `x`, or from the fitted model.

- nlag:

  Numeric vector giving the number of seasonal multiples of each period
  to test, i.e. lags `season * (1:nlag)`. A single value is recycled
  across `season`; otherwise give one value per period, e.g. `c(5, 3)`.
  Default `5`.

- lags_short:

  Number of low-order, non-seasonal lags `1:lags_short` added to the
  seasonal ones. Default `0` (none).

- npar:

  Number of estimated model parameters. Default `0` for a raw series;
  for a fitted
  [msarima](https://danigiro.github.io/msarima/reference/msarima.md)
  model it defaults to the number of estimated (non-fixed) AR, seasonal
  AR, MA and seasonal MA coefficients.

## Value

An object of class `"htest"`, with:

- statistic:

  the test statistic ("X-squared").

- parameter:

  the degrees of freedom ("df"), the number of lags tested minus `npar`.

- p.value:

  the p-value, from the upper tail of a chi-squared distribution with
  `parameter` degrees of freedom.

- method:

  a character string describing the test.

- data.name:

  a character string naming the data.

- lag:

  the lags actually tested.

## Details

The lag set tested is
`sort(unique(c(1:lags_short, outer(season, 1:nlag))))`, the same
construction as in the original two-season versions of these tests,
generalised to any number of periods.

As in [sacf](https://danigiro.github.io/msarima/reference/sacf.md),
leading and trailing `NA`s (such as those produced by differencing in
[msarima](https://danigiro.github.io/msarima/reference/msarima.md)) are
dropped
([stats::na.contiguous](https://rdrr.io/r/stats/na.contiguous.html))
before computing the statistic, and the number of observations used is
the length of the series that remains.

## References

Grigoletto, M. and Lisi, F. (2026). Portmanteau goodness-of-fit tests in
the presence of simple, mixed, and multiple seasonal autocorrelation.
*arXiv*.
[doi:10.48550/arXiv.2609.32381](https://doi.org/10.48550/arXiv.2609.32381)

Ljung, G. M. and Box, G. E. P. (1978). On a measure of lack of fit in
time series models. *Biometrika*, 65(2), 297-303.

Monti, A. C. (1994). A proposal for a residual autocorrelation test in
linear models. *Biometrika*, 81(4), 776-780.

## See also

[sacf](https://danigiro.github.io/msarima/reference/sacf.md),
[msarima](https://danigiro.github.io/msarima/reference/msarima.md),
stats::Box.test

## Examples

``` r
fit <- msarima(MAD_m40, order = c(1, 0, 2),
               seasonal = list(c(0, 1, 1), c(0, 1, 1)))

# periods 24 and 168 taken from the msts object
msportmanteau(MAD_m40, nlag = 3)
#> 
#>  Multi-seasonal Ljung-Box test
#> 
#> data:  MAD_m40
#> X-squared = 4767.5, df = 6, p-value < 2.2e-16
#> 
# residual test, periods and npar from the model
msportmanteau(fit, type = "Monti", nlag = 3)
#> 
#>  Multi-seasonal Monti test
#> 
#> data:  MAD_m40
#> X-squared = 3.8858, df = 1, p-value = 0.0487
#> 
```
