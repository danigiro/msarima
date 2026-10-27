#' Forecast from a Multi-Seasonal ARIMA Model
#'
#' @description
#' Returns forecasts and prediction intervals for a fitted multi-seasonal ARIMA model.
#'
#' @param object An object of class \code{"msarima"}.
#' @param h Number of periods to forecast. Default is 10.
#' @param level Confidence levels for the prediction intervals, in percent.
#'   Default is \code{c(80, 95)}.
#' @param fan If \code{TRUE}, \code{level} is set to \code{seq(51, 99, by = 3)},
#'   which is suitable for fan plots. Default is \code{FALSE}.
#' @param ... Additional arguments (currently not used).
#'
#' @return An object of class \code{"forecast"} containing the following elements:
#' \item{method}{The name of the forecasting method (same as \code{object$method}).}
#' \item{model}{The fitted \code{"msarima"} model object.}
#' \item{level}{The confidence levels for prediction intervals (as percentages).}
#' \item{mean}{Point forecasts as a time series (ts or msts object).}
#' \item{lower}{Lower limits for prediction intervals. Matrix with one column for each
#'   confidence level specified in \code{level}.}
#' \item{upper}{Upper limits for prediction intervals. Matrix with one column for each
#'   confidence level specified in \code{level}.}
#' \item{x}{The original time series.}
#' \item{series}{The name of the series.}
#' \item{fitted}{Fitted values from the model.}
#' \item{residuals}{Residuals from the model.}
#'
#' @details
#' Point forecasts are computed by running the model recursion forward, setting the
#' future innovations to zero. Differencing and the constant term, when present, are
#' taken from the fitted model, so the forecasts are on the scale of the original
#' series.
#'
#' Prediction intervals assume Gaussian innovations. The forecast variance
#' \eqn{V_h} is obtained from the MA(\eqn{\infty}) representation of the model.
#'
#' @examples
#' # Simulate and fit model
#' set.seed(123)
#' y <- msarima.sim(nsim = 200,
#'                  order = c(1, 0, 0),
#'                  seasonal = list(c(1, 0, 0), c(0, 0, 1)),
#'                  seasonal.period = c(4, 11),
#'                  ar = c(0.5, 0.5), ma = 0.3, sigma2 = 1)
#'
#' fit <- msarima(y, order = c(1, 0, 0),
#'                seasonal = list(c(1, 0, 0), c(0, 0, 1)),
#'                seasonal.period = c(4, 11))
#'
#' # Forecast 20 steps ahead with default confidence levels (80%, 95%)
#' fc <- forecast(fit, h = 20)
#' print(fc)
#'
#' # Forecast with custom prediction intervals
#' fc2 <- forecast(fit, h = 20, level = c(90, 95, 99))
#'
#' @seealso \code{\link{msarima}}, \code{\link{simulate.msarima}}
#'
#' @export forecast.msarima
#' @export
forecast.msarima <- function(
  object,
  h = 10,
  level = c(80, 95),
  fan = FALSE,
  ...
) {
  if (fan) {
    level <- seq(51, 99, by = 3)
  } else {
    level <- sort(level)
  }

  order_mat <- do.call(rbind, c(list(object$order), object$seasonal))
  S_vec <- c(1, object$seasonal.period) # Include 1 for non-seasonal
  n_seasonal <- length(object$seasonal)

  sigma <- sqrt(object$sigma2)

  coef_names <- names(object$coef)
  n_params <- length(object$coef)

  n_ar <- sum(order_mat[, 1])
  n_ma <- sum(order_mat[, 3])

  if (n_ar > 0) {
    ar_coefs <- object$coef[1:n_ar]
  }
  if (n_ma > 0) {
    ma_coefs <- object$coef[(n_ar + 1):(n_ar + n_ma)]
  }

  # Intercept, only present if the model was fitted with include.constant = TRUE
  if ("intercept" %in% names(object$coef)) {
    mu <- unname(object$coef[["intercept"]])
  } else {
    mu <- 0
  }

  poly_ar <- build_lag_poly(ar_coefs, order_mat[, 1], S_vec, -1)
  poly_ma <- build_lag_poly(ma_coefs, order_mat[, 3], S_vec, +1)

  if (sum(order_mat[, 2]) > 0) {
    poly_i <- do.call(
      prod,
      mapply(
        function(x, y) {
          if (y > 0) {
            polynomial(c(1, rep(0, x - 1), -1))^y
          } else {
            polynomial(1)
          }
        },
        x = S_vec,
        y = order_mat[, 2],
        SIMPLIFY = FALSE
      )
    )
  } else {
    poly_i <- polynomial(1)
  }

  poly_ar_tot <- poly_ar
  poly_ma_tot <- poly_ma
  poly <- poly_ar * poly_i

  if (length(poly_ar_tot) > 1) {
    para_ar <- -poly_ar_tot[2:length(poly_ar_tot)]
  } else {
    para_ar <- 0
  }

  if (length(poly_ma_tot) > 1) {
    para_ma <- poly_ma_tot[2:length(poly_ma_tot)]
  } else {
    para_ma <- 0
  }

  if (length(poly) > 1) {
    para <- -poly[2:length(poly)]
  } else {
    para <- 0
  }

  k_ar <- length(para_ar)
  k_ma <- length(para_ma)
  k <- length(para)

  x <- object$x
  eps <- object$residuals
  nx <- length(x)

  fcast <- numeric(h)
  for (i in 1:h) {
    ar_contrib <- if (k > 0) sum(para * x[nx:(nx - k + 1)]) else 0
    ma_contrib <- if (k_ma > 0) sum(para_ma * eps[nx:(nx - k_ma + 1)]) else 0

    fcast[i] <- mu + ar_contrib + ma_contrib

    x <- c(x, fcast[i])
    eps <- c(eps, 0)
    nx <- length(x)
  }

  # Forecast error variance for levels Y
  sigma2 <- object$sigma2

  # psi_0..psi_{h-1} for stationary ARMA part: Theta(B)/Phi(B)
  if (h == 1) {
    psi <- 1
  } else {
    psi <- c(
      1,
      stats::ARMAtoMA(
        ar = if (k_ar > 0) para_ar else numeric(0),
        ma = if (k_ma > 0) para_ma else numeric(0),
        lag.max = h - 1
      )
    )
  }

  # Delta(B) coefficients from differencing polynomial poly_i: delta_0..delta_m
  delta <- as.numeric(poly_i)
  if (abs(delta[1] - 1) > 1e-8) {
    delta <- delta / delta[1]
  }
  m <- length(delta) - 1

  # g(B) = 1/Delta(B): g_0..g_{h-1}
  g <- numeric(h)
  g[1] <- 1
  if (m > 0 && h > 1) {
    for (n in 2:h) {
      i_max <- min(m, n - 1)
      g[n] <- -sum(delta[2:(i_max + 1)] * g[(n - 1):(n - i_max)])
    }
  }

  # w(B)=g(B)*psi(B): w_0..w_{h-1}
  w_full <- convolve(g, rev(psi), type = "open")
  w <- w_full[1:h]

  # Var(e(h)) = sigma2 * sum_{n=0}^{h-1} w_n^2  (cumulative by horizon)
  fvar <- sigma2 * cumsum(w^2)

  nconf <- length(level)
  lower <- upper <- matrix(NA, nrow = h, ncol = nconf)

  for (i in 1:nconf) {
    qq <- qnorm(0.5 + level[i] / 200)
    lower[, i] <- fcast - qq * sqrt(fvar)
    upper[, i] <- fcast + qq * sqrt(fvar)
  }

  colnames(lower) <- colnames(upper) <- paste0(level, "%")

  out <- structure(
    list(
      method = object$method,
      model = object,
      level = level,
      mean = ahead_ts(object$x, fcast),
      lower = ahead_ts(object$x, lower),
      upper = ahead_ts(object$x, upper),
      x = object$x,
      series = object$series,
      fitted = object$fitted,
      residuals = object$residuals
    ),
    class = "forecast"
  )

  return(out)
}
