#' Simulate from a Multi-Seasonal ARIMA Model
#'
#' @description
#' Returns a simulated time series from a fitted ARIMA model with one or more seasonal
#' components.
#'
#' @param object An object of class \code{"msarima"}.
#' @param nsim Number of periods for the simulated series. Default is 100.
#' @param seed Either \code{NULL} or an integer that will be used in a call to
#'   \code{\link[base]{set.seed}} before simulating the time series. Default is \code{NULL}.
#' @param burn Number of burn-in observations to discard. Default is 100.
#' @param ... Additional arguments (currently not used).
#'
#' @return A time series object (ts or msts) containing the simulated series of length \code{nsim}.
#'   The time series structure matches that of the original data used to fit the model.
#'
#' @details
#' This function extracts the model specification and parameters from a fitted
#' \code{msarima} object and simulates a time series by calling \code{\link{msarima.sim}}.
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' y <- msarima.sim(nsim = 200,
#'                  order = c(1, 0, 0),
#'                  seasonal = list(c(1, 0, 0), c(0, 0, 1)),
#'                  seasonal.period = c(4, 11),
#'                  ar = c(0.5, 0.5), ma = 0.3, sigma2 = 1)
#'
#' # Fit model
#' fit <- msarima(y, order = c(1, 0, 0),
#'                seasonal = list(c(1, 0, 0), c(0, 0, 1)),
#'                seasonal.period = c(4, 11))
#'
#' # Simulate from fitted model
#' y_sim <- simulate(fit, nsim = 100)
#'
#' @seealso \code{\link{msarima}}, \code{\link{msarima.sim}}
#' @method simulate msarima
#'
#' @export
simulate.msarima <- function(
  object = NULL,
  nsim = 100,
  seed = NULL,
  burn = 100,
  ...
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Extract model specification
  order_mat <- do.call(rbind, c(list(object$order), object$seasonal))
  sum_order <- colSums(order_mat)[-2]

  # Extract coefficients
  ar_idx <- if (sum_order[1] != 0) 1:sum_order[1] else NULL
  ma_idx <- if (sum_order[2] != 0) sum_order[1] + (1:sum_order[2]) else NULL

  sigma2 <- object$sigma2

  # Call msarima.sim with extracted parameters
  msarima.sim(
    nsim = nsim,
    order = object$order,
    seasonal = object$seasonal,
    seasonal.period = object$seasonal.period,
    ar = if (length(ar_idx) > 0) object$coef[ar_idx] else NULL,
    ma = if (length(ma_idx) > 0) object$coef[ma_idx] else NULL,
    sigma2 = sigma2,
    burn = burn
  )
}

#' Simulate from a Multi-Seasonal ARIMA Model (Direct Interface)
#'
#' @description
#' Simulates a time series from a multi-seasonal ARIMA model with specified parameters.
#' This function provides a direct interface for simulation without requiring a fitted model object.
#'
#' @param nsim Number of periods for the simulated series (required).
#' @param order A numeric vector of length 3 specifying the non-seasonal part of the
#'   ARIMA model: the AR order (p), the degree of differencing (d), and the MA order (q).
#'   Default is \code{c(0, 0, 0)}.
#' @param seasonal A list of numeric vectors, each of length 3, specifying the seasonal
#'   components. Each vector gives the seasonal AR order (P), seasonal differencing order (D),
#'   and seasonal MA order (Q) for that seasonal component.
#'   The list must contain at least one seasonal component.
#'   Default is \code{list(c(0, 0, 0))}.
#' @param seasonal.period A numeric vector specifying the seasonal periods corresponding
#'   to each element in \code{seasonal}. Default is \code{NULL}.
#' @param ar Numeric vector of AR coefficients (non-seasonal and seasonal combined).
#'   Length must match the total AR order from \code{order} and all \code{seasonal} components.
#'   Default is \code{NULL}.
#' @param ma Numeric vector of MA coefficients (non-seasonal and seasonal combined).
#'   Length must match the total MA order from \code{order} and all \code{seasonal} components.
#'   Default is \code{NULL}.
#' @param sigma2 Innovation variance. Default is 1.
#' @param burn Number of burn-in observations to discard. If \code{NULL}, automatically
#'   set to \code{max(100, 2 * kk)} where \code{kk} is the maximum lag in the model.
#'   Default is \code{NULL}.
#'
#' @return A time series object (msts) containing the simulated series of length \code{nsim}.
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' y <- msarima.sim(nsim = 200,
#'                  order = c(1, 0, 0),
#'                  seasonal = list(c(1, 0, 0), c(0, 0, 1)),
#'                  seasonal.period = c(4, 11),
#'                  ar = c(0.5, 0.5), ma = 0.3, sigma2 = 1)
#'
#' @seealso \code{\link{msarima}}, \code{\link{simulate.msarima}}
#' @export
msarima.sim <- function(
  nsim,
  order = c(0, 0, 0),
  seasonal = list(c(0, 0, 0), c(0, 0, 0)),
  seasonal.period = NULL,
  ar = NULL,
  ma = NULL,
  sigma2 = 1,
  burn = NULL
) {
  # Validate seasonal input
  if (!is.list(seasonal)) {
    stop("'seasonal' must be a list")
  }
  n_seasonal <- length(seasonal)

  if (n_seasonal < 2) {
    stop("'seasonal' must have at least two seasonal components")
  }

  # Check that each seasonal element is a vector of length 3
  for (i in 1:n_seasonal) {
    if (!is.numeric(seasonal[[i]]) || length(seasonal[[i]]) != 3) {
      stop(sprintf("seasonal[[%d]] must be a numeric vector of length 3", i))
    }
    if (any(seasonal[[i]] < 0) || any(seasonal[[i]] != floor(seasonal[[i]]))) {
      stop(sprintf("seasonal[[%d]] must contain non-negative integers", i))
    }
  }

  # Build order matrix (same structure as in msarima)
  order_mat <- do.call(rbind, c(list(order), seasonal))
  S_vec <- c(1, seasonal.period) # Include 1 for non-seasonal

  # Count total AR and MA parameters
  n_ar <- sum(order_mat[, 1])
  n_ma <- sum(order_mat[, 3])

  # Validate parameters match orders
  if (n_ar > 0 && (is.null(ar) || length(ar) != n_ar)) {
    stop(sprintf("Length of AR coefficients must be %d to match orders", n_ar))
  }
  if (n_ma > 0 && (is.null(ma) || length(ma) != n_ma)) {
    stop(sprintf("Length of MA coefficients must be %d to match orders", n_ma))
  }

  sigma <- sqrt(sigma2)

  # Build AR polynomial
  if (n_ar > 0) {
    ar_comp <- rep(S_vec, order_mat[, 1])
    ar_beta <- tapply(ar, ar_comp, function(x) x, simplify = FALSE)

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
    ma_beta <- tapply(ma, ma_comp, function(x) x, simplify = FALSE)

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
  kk <- max(k_ar, k_ma, k)

  # Determine burn-in
  if (is.null(burn)) {
    burn <- max(100, 2 * kk)
  }

  n_total <- nsim + burn
  # Generate innovations
  y <- rnorm(n_total, mean = 0, sd = sigma)
  eps <- rnorm(n_total, mean = 0, sd = sigma)

  # Generate the series (exact same loop as original msarima.sim)
  for (j in (kk + 1):n_total) {
    y[j] <- sum(para * y[(j - 1):(j - k)]) +
      sum(para_ma * eps[(j - 1):(j - k_ma)]) +
      eps[j]
  }

  # Remove burn-in
  y_sim <- y[(burn + 1):n_total]

  # Create ts or msts object
  valid_periods <- seasonal.period[!is.na(seasonal.period)]
  if (length(valid_periods) > 1) {
    y_sim <- msts(y_sim, seasonal.periods = valid_periods)
  } else if (length(valid_periods) == 1) {
    y_sim <- ts(y_sim, frequency = valid_periods[1])
  } else {
    y_sim <- ts(y_sim)
  }

  return(y_sim)
}
