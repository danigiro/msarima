# Getting started with msarima

Hourly, daily or half-hourly data usually show more than one seasonal
cycle at the same time: the number of cars on a motorway repeats itself
every 24 hours *and* every 168 hours (a week).
[`stats::arima()`](https://rdrr.io/r/stats/arima.html) can accommodate
only one of the two. The `msarima` package fits the multiple seasonal
model of Lisi and Grigoletto (2025),

\\ \phi(B)\\\Phi_1(B^{S_1})\cdots\Phi_k(B^{S_k})\\ \nabla^d
\nabla\_{S_1}^{D_1}\cdots\nabla\_{S_k}^{D_k}\\ y_t =
\theta(B)\\\Theta_1(B^{S_1})\cdots\Theta_k(B^{S_k})\\ e_t , \\

with one autoregressive, differencing and moving average operator for
each seasonal period \\S_1, \dots, S_k\\.

``` r

library(msarima)
```

## The data

`MAD_m40` is the hourly number of vehicles recorded by a sensor on the
M-40 ring road of Madrid between February and May 2024. It is an `msts`
object with seasonal periods 24 and 168 (2904 observations in all), so
those periods do not have to be repeated in every call.

``` r

head(MAD_m40, 12)
#> Multi-Seasonal Time Series:
#> Start: 1 1
#> Seasonal Periods: 24 168
#> Data:
#>  [1]  363  214  175  116  112  103  179  514 1241 1306 1127 1053

plot(window(MAD_m40, end = 4), ylab = "vehicles per hour",
     main = "M-40, first four weeks (time in weeks)")
```

![Line chart of vehicles per hour over the first four weeks. A daily
cycle of two peaks and a night trough repeats, with lower peaks on the
last two days of every week.](msarima_files/figure-html/data-1.png)

Both cycles are visible: the daily peak-hour pattern, and the weekly one
with the weekends clearly lower.

## Looking at the cycles

[`sacf()`](https://danigiro.github.io/msarima/reference/sacf.md) and
[`spacf()`](https://danigiro.github.io/msarima/reference/sacf.md) are
the correlograms of the package. They work like
[`stats::acf()`](https://rdrr.io/r/stats/acf.html), but the lags
belonging to each seasonal period are coloured, so the cycles can be
read off the plot instead of counted by hand. The periods are taken from
the `msts` object.

``` r

sacf(MAD_m40, lag.max = 200)
```

![Correlogram up to lag 200. The autocorrelations stay well above the
confidence band and peak at the multiples of 24, in orange, and at 168,
in blue, decaying very slowly.](msarima_files/figure-html/acf-1.png)

``` r

spacf(MAD_m40, lag.max = 200)
```

![Partial correlogram up to lag 200. Only the first lags and those
around the multiples of 24 and 168 stand out of the confidence
band.](msarima_files/figure-html/pacf-1.png)

Orange marks the multiples of 24, blue the multiples of 168; `col`
changes them. The autocorrelation decays very slowly at both, which
suggests one difference at each period.

There is a `ggplot2` version of both,
[`ggsacf()`](https://danigiro.github.io/msarima/reference/sacf.md) and
[`ggspacf()`](https://danigiro.github.io/msarima/reference/sacf.md):

``` r

ggsacf(MAD_m40, lag.max = 200)
```

![The same correlogram in ggplot2, with a legend mapping orange to
period 24 and blue to period
168.](msarima_files/figure-html/ggacf-1.png)

## Fitting the model

The non-seasonal part goes in `order`, and each seasonal component is
one element of the `seasonal` list, in the same order as the periods:

``` r

fit <- msarima(MAD_m40,
               order = c(1, 0, 2),
               seasonal = list(c(0, 1, 1), c(0, 1, 1)))
fit
#> Series: MAD_m40 
#> MSARIMA(1,0,2)(0,1,1)[24](0,1,1)[168]
#> 
#> Coefficients:
#>    AR(1)    MA(1)    MA(2) SMA(1,1) SMA(2,1) 
#>   0.8235  -0.0849  -0.0629  -0.7518  -0.8184 
#> 
#> sigma^2 = 5105.4386  log likelihood = -12933.63
#> AIC = 25879.25   AICc = 25879.28   BIC = 25914.69
```

[`summary()`](https://rdrr.io/r/base/summary.html) adds standard errors,
obtained from the numerical Hessian of the log-likelihood, and the
associated tests:

``` r

summary(fit)
#> Series: MAD_m40 
#> MSARIMA(1,0,2)(0,1,1)[24](0,1,1)[168]
#> 
#> Coefficients:
#>          Estimate Std. Error z value Pr(>|z|)    
#> AR(1)     0.82355    0.02055  40.082  < 2e-16 ***
#> MA(1)    -0.08489    0.03043  -2.790  0.00527 ** 
#> MA(2)    -0.06288    0.02893  -2.173  0.02976 *  
#> SMA(1,1) -0.75184    0.01762 -42.681  < 2e-16 ***
#> SMA(2,1) -0.81843    0.01107 -73.901  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> sigma^2 = 5105.4386  log likelihood = -12933.6262
#> AIC = 25879.2524   AICc = 25879.2834   BIC = 25914.6850
```

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) are all available and behave
as usual. When the model contains differencing the series is differenced
before estimation, so the first \\d + \sum_k D_k S_k\\ residuals are
`NA`:

``` r

c(nobs = fit$nobs, n.used = fit$n.used, leading_NA = sum(is.na(residuals(fit))))
#>       nobs     n.used leading_NA 
#>       2904       2712        192
```

Parameters can be held fixed by supplying `fixed` with one entry per
parameter and `NA` where the parameter is to be estimated. The order is
AR, seasonal AR, MA, seasonal MA, sigma:

``` r

fit_fixed <- msarima(MAD_m40,
                     order = c(1, 0, 2),
                     seasonal = list(c(0, 1, 1), c(0, 1, 1)),
                     fixed = c(NA, NA, NA, -0.6, NA, NA))
c(free = AIC(fit), fixed = AIC(fit_fixed))
#>     free    fixed 
#> 25879.25 25930.09
```

## Checking the residuals

Passing the fitted model to
[`sacf()`](https://danigiro.github.io/msarima/reference/sacf.md) uses
the residuals and the seasonal periods of the model, so no argument is
needed. What is left at the seasonal lags is what the model has failed
to capture:

``` r

sacf(fit, lag.max = 200, main = "Residuals")
```

![Correlogram of the residuals up to lag 200. The spikes are much
smaller than before and mostly inside the confidence
band.](msarima_files/figure-html/resid-1.png)

## Forecasting

[`forecast()`](https://generics.r-lib.org/reference/forecast.html)
returns an object of the usual class of the **forecast** package, so it
prints and plots like any other. Prediction intervals come from the
\\\mathrm{MA}(\infty)\\ representation of the model, under normal
innovations:

``` r

fc <- forecast(fit, h = 48)
plot(fc, include = 336, main = "Two days ahead", ylab = "vehicles per hour")
```

![Line chart of the last two weeks of data followed by 48 hours of
forecasts, which repeat the daily shape, inside shaded 80 and 95 per
cent prediction bands.](msarima_files/figure-html/forecast-1.png)

## Simulating

[`msarima.sim()`](https://danigiro.github.io/msarima/reference/msarima.sim.md)
generates a series from a model specification, without any data:

``` r

set.seed(1)
y <- msarima.sim(nsim = 500,
                 order = c(1, 0, 1),
                 seasonal = list(c(0, 0, 1), c(0, 0, 1)),
                 seasonal.period = c(7, 30),
                 ar = 0.6, ma = c(0.4, 0.3, 0.2),
                 sigma2 = 0.5)
```

and [`simulate()`](https://rdrr.io/r/stats/simulate.html) does the same
from a fitted model, reusing its specification and its estimates:

``` r

y_sim <- simulate(fit, nsim = 168, seed = 1)
```

## References

Lisi, F. and Grigoletto, M. (2025). Multiple seasonal autoregressive
integrated moving average models. *Journal of Forecasting*, **44**(6),
2037–2052.
