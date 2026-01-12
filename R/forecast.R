#' Forecast from Multi-Seasonal ARIMA Models
#'
#' @description
#' Returns forecasts and other information for a fitted multi-seasonal ARIMA model.
#'
#' @param object An object of class \code{"msarima"}.
#' @param h Number of periods for forecasting. Default is 10.
#' @param level Confidence levels for prediction intervals.
#' @param fan If \code{TRUE}, \code{level} is set to \code{seq(51, 99, by = 3)}.
#'   This is suitable for fan plots.
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
#' The forecast is computed recursively using the fitted MSARIMA model.
#' Prediction intervals are computed assuming normally distributed errors and
#' the forecast variance \eqn{V_h} is computed using the MA(\eqn{\infty}) representation
#' of the model.
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
  # Handle fan parameter
  if (fan) {
    level <- seq(51, 99, by = 3)
  } else {
    level <- sort(level)
  }

  # Extract model specification as matrix (same structure as in msarima)
  order_mat <- do.call(rbind, c(list(object$order), object$seasonal))
  S_vec <- c(1, object$seasonal.period) # Include 1 for non-seasonal
  n_seasonal <- length(object$seasonal)

  sigma <- sqrt(object$sigma2)

  # Extract coefficients
  coef_names <- names(object$coef)
  n_params <- length(object$coef)

  # Count AR and MA parameters
  n_ar <- sum(order_mat[, 1])
  n_ma <- sum(order_mat[, 3])

  # Split coefficients into AR and MA
  if (n_ar > 0) {
    ar_coefs <- object$coef[1:n_ar]
  }
  if (n_ma > 0) {
    ma_coefs <- object$coef[(n_ar + 1):(n_ar + n_ma)]
  }

  # Build AR polynomial
  if (n_ar > 0) {
    ar_comp <- rep(S_vec, order_mat[, 1])
    ar_beta <- tapply(ar_coefs, ar_comp, function(x) x, simplify = FALSE)

    poly_ar <- do.call(
      prod,
      mapply(
        function(x, y) polynomial(c(1, rep(0, x - 1), -y)),
        x = unique(ar_comp),
        y = ar_beta,
        SIMPLIFY = FALSE
      )
    )
  } else {
    poly_ar <- polynomial(1)
  }

  # Build MA polynomial
  if (n_ma > 0) {
    ma_comp <- rep(S_vec, order_mat[, 3])
    ma_beta <- tapply(ma_coefs, ma_comp, function(x) x, simplify = FALSE)

    poly_ma <- do.call(
      prod,
      mapply(
        function(x, y) polynomial(c(1, rep(0, x - 1), y)),
        x = unique(ma_comp),
        y = ma_beta,
        SIMPLIFY = FALSE
      )
    )
  } else {
    poly_ma <- polynomial(1)
  }

  # Build integration polynomial
  if (sum(order_mat[, 2]) > 0) {
    poly_i <- do.call(
      prod,
      mapply(
        function(x, y) {
          if (y > 0) {
            polynomial(c(1, rep(0, x - 1), -1))
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

  # Total polynomials
  poly_ar_tot <- poly_ar
  poly_ma_tot <- poly_ma
  poly <- poly_ar * poly_i

  # Extract parameters
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

  # Get data and residuals
  x <- object$x
  eps <- object$residuals
  nx <- length(x)

  # Compute point forecasts
  fcast <- numeric(h)
  for (i in 1:h) {
    ar_contrib <- if (k > 0) sum(para * x[nx:(nx - k + 1)]) else 0
    ma_contrib <- if (k_ma > 0) sum(para_ma * eps[nx:(nx - k_ma + 1)]) else 0

    fcast[i] <- ar_contrib + ma_contrib

    # Update series and errors for next forecast
    x <- c(x, fcast[i])
    eps <- c(eps, 0)
    nx <- length(x)
  }

  # ---- Forecast error variance for levels Y ----
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
  if (m > 0) {
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

  # Compute prediction intervals
  nconf <- length(level)
  lower <- upper <- matrix(NA, nrow = h, ncol = nconf)

  for (i in 1:nconf) {
    qq <- qnorm(0.5 + level[i] / 200)
    lower[, i] <- fcast - qq * sqrt(fvar)
    upper[, i] <- fcast + qq * sqrt(fvar)
  }

  colnames(lower) <- colnames(upper) <- paste0(level, "%")

  # Return forecast object
  out <- structure(
    list(
      method = object$method,
      model = object,
      level = level,
      mean = future_msts(object$x, fcast),
      lower = future_msts(object$x, lower),
      upper = future_msts(object$x, upper),
      x = object$x,
      series = object$series,
      fitted = object$fitted,
      residuals = object$residuals
    ),
    class = "forecast"
  )

  return(out)
}

# Copy msts attributes from x to y shifted to forecast period
# Version from the forecast package:
# https://github.com/robjhyndman/forecast/blob/29754e2e8044c7b357a0cae321ff34df24ec9386/R/msts.R#L116
future_msts <- function(x, y) {
  if (NCOL(y) > 1) {
    class(y) <- c("mts", "ts", "matrix")
  } else {
    class(y) <- "ts"
  }
  if (inherits(x, "msts")) {
    class(y) <- c("msts", class(y))
  }
  attr <- attributes(x)
  attr$tsp[1:2] <- attr$tsp[2] + c(1, NROW(y)) / attr$tsp[3]
  attributes(y)$tsp <- attr$tsp
  attributes(y)$msts <- attr$msts
  y
}

# Copy msts attributes from x to y
# Version from the forecast package:
# https://github.com/robjhyndman/forecast/blob/29754e2e8044c7b357a0cae321ff34df24ec9386/R/msts.R#L90
copy_msts <- function(x, y) {
  if (NROW(x) > NROW(y)) {
    # Pad y with initial NAs
    if (NCOL(y) == 1) {
      y <- c(rep(NA, NROW(x) - NROW(y)), y)
    } else {
      y <- rbind(matrix(NA, ncol = NCOL(y), nrow = NROW(x) - NROW(y)), y)
    }
  } else if (NROW(x) != NROW(y)) {
    stop("x and y should have the same number of observations")
  }
  if (NCOL(y) > 1) {
    class(y) <- c("mts", "ts", "matrix")
  } else {
    class(y) <- "ts"
  }
  if (inherits(x, "msts")) {
    class(y) <- c("msts", class(y))
  }
  attr <- attributes(x)
  attributes(y)$tsp <- attr$tsp
  attributes(y)$msts <- attr$msts
  y
}
