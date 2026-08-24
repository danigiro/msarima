# Fit a Multi-Seasonal ARIMA Model

Fits an ARIMA model with one or more seasonal components to a univariate
time series. This function extends the standard ARIMA framework to time
series with several seasonal patterns, such as daily and weekly
seasonality in hourly data, or weekly and yearly seasonality in daily
data.

## Usage

``` r
msarima(
  y,
  order = c(0, 0, 0),
  seasonal = list(),
  seasonal.period = NULL,
  fixed = NULL,
  include.constant = FALSE,
  optim_control = list(maxit = 200)
)
```

## Arguments

- y:

  A numeric vector or time series object (`ts` or `msts`) containing the
  data to be modelled. If `y` is a `ts` object and `seasonal.period` is
  not supplied, its frequency is used as the first seasonal period; if
  `y` is an `msts` object, its seasonal periods are used.

- order:

  A numeric vector of length 3 specifying the non-seasonal part of the
  ARIMA model: the AR order (p), the degree of differencing (d), and the
  MA order (q). Default is `c(0, 0, 0)`.

- seasonal:

  A list of numeric vectors, each of length 3, specifying the seasonal
  components. Each vector gives the seasonal AR order (P), the seasonal
  differencing order (D), and the seasonal MA order (Q) of one
  component. For example, `list(c(1, 0, 1), c(0, 1, 1))` specifies two
  seasonal components. Use [`list()`](https://rdrr.io/r/base/list.html)
  for a model without seasonal components. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- seasonal.period:

  A numeric vector of seasonal periods, one for each element of
  `seasonal`, for example `c(7, 30)` for weekly and monthly seasonality.
  If `NULL`, the periods are taken from `y` when it is a `ts` or `msts`
  object. Default is `NULL`.

- fixed:

  A numeric vector with one entry per parameter, used to hold some of
  them at a given value: only the `NA` entries are estimated. The
  parameters are ordered as non-seasonal AR, seasonal AR (one block per
  seasonal component), non-seasonal MA, seasonal MA (one block per
  seasonal component), sigma, and the constant term when
  `include.constant = TRUE`. Default is `NULL`, meaning that all
  parameters are estimated.

- include.constant:

  Should the model include a constant term? The constant is the
  intercept of the model equation, not the mean of the series: for a
  stationary model the two are related by
  `mean = intercept / prod(1 - ar)`. A constant cannot be combined with
  differencing and is dropped, with a warning, when `d` or any `D` is
  non-zero. Default is `FALSE`.

- optim_control:

  A list of control parameters passed to
  [`optim`](https://rdrr.io/r/stats/optim.html). Default is
  `list(maxit = 200)`.

## Value

An object of class `"msarima"`, which is a list containing:

- order:

  The non-seasonal ARIMA order (p, d, q).

- seasonal:

  List of seasonal component specifications.

- seasonal.period:

  Vector of seasonal periods.

- coef:

  Named vector of estimated coefficients: the AR, SAR, MA and SMA
  parameters, followed by the intercept if `include.constant = TRUE`.

- var_coef:

  Variance-covariance matrix of the coefficient estimates.

- mask:

  Logical vector indicating which parameters were fixed.

- loglik:

  Log-likelihood of the fitted model.

- aic:

  Akaike Information Criterion.

- bic:

  Bayesian Information Criterion.

- aicc:

  Corrected Akaike Information Criterion.

- sigma2:

  Estimated innovation variance.

- x:

  The original time series.

- fitted:

  One-step-ahead fitted values, on the scale of the original series.

- residuals:

  Model residuals. Together with `fitted`, the first `d + sum(D * S)`
  values are `NA`, since differencing leaves no observation to predict
  there.

- method:

  Description of the fitted model.

- series:

  Name of the series.

- nobs:

  Number of observations in the original series.

- n.used:

  Number of observations used in estimation (after differencing).

- call:

  The matched call.

- convergence:

  Convergence code from `optim` (0 indicates successful convergence).

## Details

The model is estimated by conditional maximum likelihood: the residuals
are computed recursively and the resulting Gaussian log-likelihood is
maximised numerically with
[`optim`](https://rdrr.io/r/stats/optim.html), using the BFGS algorithm.

When `d` or any `D` is non-zero, the series is differenced before
estimation and the ARMA part is estimated on the differenced data. The
first `d + sum(D * S)` observations are therefore not used, which is why
the residuals start with as many `NA` values.

Standard errors are computed from the inverse of the observed
information matrix (negative Hessian). If the Hessian is singular, a
warning is issued and the standard errors are set to `NA`.

## Note

- The series must have at least 10 observations, both before and after
  differencing.

- Seasonal periods must be at least 2 and at most half the length of the
  series.

- If `y` is a `ts` or `msts` object, `x`, `fitted` and `residuals` keep
  the same time series structure. If `y` is a plain numeric vector, they
  are returned as `ts` with one seasonal period, or as `msts` with
  several.

## Examples

``` r
# Example 1: ARIMA with two seasonal components
set.seed(123)
y <- msarima.sim(nsim = 200,
                 order = c(1, 0, 0),
                 seasonal = list(c(1, 0, 0), c(0, 0, 1)),
                 seasonal.period = c(4, 11),
                 ar = c(0.5, 0.5), ma = 0.3, sigma2 = 1)

fit <- msarima(y, order = c(1, 0, 0),
               seasonal = list(c(1, 0, 0), c(0, 0, 1)),
               seasonal.period = c(4, 11))
print(fit)
#> Series: y 
#> MSARIMA(1,0,0)(1,0,0)[4](0,0,1)[11]
#> 
#> Coefficients:
#>    AR(1) SAR(1,1) SMA(2,1) 
#>   0.4101   0.4678   0.3805 
#> 
#> sigma^2 = 0.9272  log likelihood = -92.44
#> AIC = 192.88   AICc = 193.08   BIC = 206.07
summary(fit)
#> Series: y 
#> MSARIMA(1,0,0)(1,0,0)[4](0,0,1)[11]
#> 
#> Coefficients:
#>          Estimate Std. Error z value Pr(>|z|)    
#> AR(1)     0.41006    0.06570   6.241 4.35e-10 ***
#> SAR(1,1)  0.46782    0.06428   7.278 3.40e-13 ***
#> SMA(2,1)  0.38051    0.06943   5.481 4.24e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> sigma^2 = 0.9272  log likelihood = -92.4387
#> AIC = 192.8775   AICc = 193.0826   BIC = 206.0707

# Example 2: msts input, seasonal periods detected automatically
y_msts <- msts(y, seasonal.periods = c(4, 11))
fit_msts <- msarima(y_msts, order = c(1, 0, 0),
                    seasonal = list(c(1, 0, 0), c(0, 0, 1)))

# Example 3: differencing
fit_diff <- msarima(cumsum(y), order = c(1, 1, 0),
                    seasonal = list(c(0, 0, 1)),
                    seasonal.period = 4)
head(residuals(fit_diff))  # the first residual is NA
#>         Qtr1       Qtr2       Qtr3       Qtr4
#> 1         NA -0.7053174 -2.4672228 -0.1094124
#> 2  0.4100991 -2.0204887                      

# Example 4: fixing parameters
# Note: the length of 'fixed' must match AR + SAR + MA + SMA + sigma
fit_fixed <- msarima(y, order = c(1, 0, 0),
                     seasonal = list(c(1, 0, 0), c(0, 0, 1)),
                     seasonal.period = c(4, 11),
                     fixed = c(0.5, NA, NA, NA))
AIC(fit); AIC(fit_fixed)
#> [1] 192.8775
#> [1] 192.7455
```
