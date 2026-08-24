# Forecast from a Multi-Seasonal ARIMA Model

Returns forecasts and prediction intervals for a fitted multi-seasonal
ARIMA model.

## Usage

``` r
# S3 method for class 'msarima'
forecast(object, h = 10, level = c(80, 95), fan = FALSE, ...)
```

## Arguments

- object:

  An object of class `"msarima"`.

- h:

  Number of periods to forecast. Default is 10.

- level:

  Confidence levels for the prediction intervals, in percent. Default is
  `c(80, 95)`.

- fan:

  If `TRUE`, `level` is set to `seq(51, 99, by = 3)`, which is suitable
  for fan plots. Default is `FALSE`.

- ...:

  Additional arguments (currently not used).

## Value

An object of class `"forecast"` containing the following elements:

- method:

  The name of the forecasting method (same as `object$method`).

- model:

  The fitted `"msarima"` model object.

- level:

  The confidence levels for prediction intervals (as percentages).

- mean:

  Point forecasts as a time series (ts or msts object).

- lower:

  Lower limits for prediction intervals. Matrix with one column for each
  confidence level specified in `level`.

- upper:

  Upper limits for prediction intervals. Matrix with one column for each
  confidence level specified in `level`.

- x:

  The original time series.

- series:

  The name of the series.

- fitted:

  Fitted values from the model.

- residuals:

  Residuals from the model.

## Details

Point forecasts are computed by running the model recursion forward,
setting the future innovations to zero. Differencing and the constant
term, when present, are taken from the fitted model, so the forecasts
are on the scale of the original series.

Prediction intervals assume normally distributed innovations. The
forecast variance \\V_h\\ is obtained from the MA(\\\infty\\)
representation of the model.

## See also

[`msarima`](https://danigiro.github.io/msarima/reference/msarima.md),
[`simulate.msarima`](https://danigiro.github.io/msarima/reference/simulate.msarima.md)

## Examples

``` r
# Simulate and fit model
set.seed(123)
y <- msarima.sim(nsim = 200,
                 order = c(1, 0, 0),
                 seasonal = list(c(1, 0, 0), c(0, 0, 1)),
                 seasonal.period = c(4, 11),
                 ar = c(0.5, 0.5), ma = 0.3, sigma2 = 1)

fit <- msarima(y, order = c(1, 0, 0),
               seasonal = list(c(1, 0, 0), c(0, 0, 1)),
               seasonal.period = c(4, 11))

# Forecast 20 steps ahead with default confidence levels (80%, 95%)
fc <- forecast(fit, h = 20)
print(fc)
#>          Point Forecast      Lo 80      Hi 80     Lo 95      Hi 95
#> 19.18182    -0.08424617 -1.3182438  1.1497514 -1.971483  1.8029903
#> 19.27273    -2.41349542 -3.7472102 -1.0797807 -4.453236 -0.3737547
#> 19.36364    -0.80138911 -2.1511477  0.5483694 -2.865667  1.2628885
#> 19.45455     0.55732814 -0.7951094  1.9097657 -1.511047  2.6257029
#> 19.54545    -0.06586351 -1.5503983  1.4186713 -2.336264  2.2045365
#> 19.63636    -1.01644566 -2.5220543  0.4891630 -3.319075  1.2861841
#> 19.72727    -0.71668606 -2.2258093  0.7924372 -3.024691  1.5913188
#> 19.81818     0.04960625 -1.4601072  1.5593197 -2.259301  2.3585137
#> 19.90909    -0.87664376 -2.4134641  0.6601766 -3.227008  1.4737201
#> 20.00000    -0.42841480 -1.9697463  1.1129167 -2.785678  1.9288482
#> 20.09091     0.36246716 -1.1796215  1.9045559 -1.995954  2.7208883
#> 20.18182     0.30932143 -1.3085505  1.9271934 -2.165000  2.7836431
#> 20.27273    -0.29278496 -1.9433733  1.3578034 -2.817142  2.2315720
#> 20.36364    -0.15231023 -1.8083362  1.5037158 -2.684983  2.3803630
#> 20.45455     0.18929555 -1.4676430  1.8462341 -2.344773  2.7233644
#> 20.54545     0.15279500 -1.5217527  1.8273427 -2.408205  2.7137947
#> 20.63636    -0.13365247 -1.8156773  1.5483724 -2.706088  2.4387826
#> 20.72727    -0.06989301 -1.7531719  1.6133859 -2.644246  2.5044599
#> 20.81818     0.08911330 -1.5943763  1.7726029 -2.485562  2.6637885
#> 20.90909     0.07170871 -1.6156150  1.7590324 -2.508830  2.6522477

# Forecast with custom prediction intervals
fc2 <- forecast(fit, h = 20, level = c(90, 95, 99))
```
