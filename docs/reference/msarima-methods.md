# Methods for msarima Objects

Methods for extracting information from and displaying fitted
multi-seasonal ARIMA models.

## Usage

``` r
# S3 method for class 'msarima'
print(x, ...)

# S3 method for class 'msarima'
summary(object, ...)

# S3 method for class 'summary_msarima'
print(x, digits = 4, signif.stars = TRUE, ...)

# S3 method for class 'msarima'
coef(object, ...)

# S3 method for class 'msarima'
fitted(object, ...)

# S3 method for class 'msarima'
residuals(object, ...)

# S3 method for class 'msarima'
logLik(object, ...)

# S3 method for class 'msarima'
AIC(object, ...)

# S3 method for class 'msarima'
BIC(object, ...)
```

## Arguments

- x, object:

  An object of class `"msarima"`.

- ...:

  Additional arguments (currently not used).

- digits:

  Number of digits to display. Default is 4.

- signif.stars:

  Logical; if `TRUE`, significance stars are printed. Default is `TRUE`.

## Value

- `print`:

  Returns the input object invisibly.

- `summary`:

  Returns an object of class `"summary_msarima"` containing:

  - `call`: The matched call.

  - `series`: Name of the time series.

  - `order`: Non-seasonal ARIMA order.

  - `seasonal`: List of seasonal component specifications.

  - `seasonal.period`: Vector of seasonal periods.

  - `coefficients`: Matrix with estimates, standard errors, z-values,
    and p-values.

  - `sigma2`: Estimated innovation variance.

  - `loglik`: Log-likelihood of the fitted model.

  - `aic`: Akaike Information Criterion.

  - `bic`: Bayesian Information Criterion.

  - `aicc`: Corrected Akaike Information Criterion.

  - `nobs`: Number of observations in the original series.

  - `n.used`: Number of observations used in estimation.

  - `method`: Description of the fitted model.

  - `convergence`: Convergence code from optimization.

- `coef`:

  A named numeric vector of coefficients.

- `fitted`:

  A time series object of fitted values.

- `residuals`:

  A time series object of residuals.

- `logLik`:

  An object of class `"logLik"`.

- `AIC`:

  The AIC value.

- `BIC`:

  The BIC value.

## Details

- `print`:

  Prints a concise summary of the fitted model including coefficients
  and information criteria.

- `summary`:

  Produces a detailed summary including coefficient estimates, standard
  errors, z-values, and p-values.

- `coef`:

  Extracts the estimated model coefficients.

- `fitted`:

  Extracts the fitted values.

- `residuals`:

  Extracts the model residuals.

- `logLik`:

  Extracts the log-likelihood value.

- `AIC`:

  Calculates the Akaike Information Criterion.

- `BIC`:

  Calculates the Bayesian Information Criterion.

## See also

[`msarima`](https://danigiro.github.io/msarima/reference/msarima.md)

## Examples

``` r
# Fit a model
set.seed(123)
y <- msarima.sim(nsim = 200,
                 order = c(1, 0, 0),
                 seasonal = list(c(1, 0, 0), c(0, 0, 1)),
                 seasonal.period = c(4, 11),
                 ar = c(0.5, 0.5), ma = 0.3, sigma2 = 1)

fit <- msarima(y, order = c(1, 0, 0),
               seasonal = list(c(1, 0, 0), c(0, 0, 1)),
               seasonal.period = c(4, 11))

# Use methods
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
coef(fit)
#>     AR(1)  SAR(1,1)  SMA(2,1) 
#> 0.4100567 0.4678163 0.3805071 
head(fitted(fit))
#> Multi-Seasonal Time Series:
#> Start: 1 1
#> Seasonal Periods: 4 11
#> Data:
#> [1] -0.8118898 -0.6042914 -0.8312793 -1.5763490 -0.7240990 -0.5768790
#> 
head(residuals(fit))
#> Multi-Seasonal Time Series:
#> Start: 1 1
#> Seasonal Periods: 4 11
#> Data:
#> [1]  0.13471646 -0.62643432 -2.21643836  0.02942345  0.40456177 -1.85157381
#> 
logLik(fit)
#> 'log Lik.' -92.43873 (df=4)
AIC(fit)
#> [1] 192.8775
BIC(fit)
#> [1] 206.0707
```
