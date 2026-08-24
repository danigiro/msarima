#' Simulate from a Multi-Seasonal ARIMA Model
#'
#' @description
#' Returns a simulated time series from a fitted ARIMA model with one or more seasonal
#' components.
#'
#' @param object An object of class \code{"msarima"}.
#' @param nsim Length of the simulated series. Default is 100.
#' @param seed Either \code{NULL} or an integer passed to \code{\link[base]{set.seed}}
#'   before simulating, to make the result reproducible. Default is \code{NULL}.
#' @param burn Number of burn-in observations to discard. Default is 100.
#' @param ... Additional arguments (currently not used).
#'
#' @return A time series object (\code{ts} or \code{msts}) of length \code{nsim}, with the
#'   seasonal periods of the fitted model.
#'
#' @details
#' The model specification and the estimated parameters, including the constant term when
#' the model was fitted with \code{include.constant = TRUE}, are extracted from the fitted
#' object and passed to \code{\link{msarima.sim}}.
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
  if (!is.null(seed)) {
    set.seed(seed)
  }

  order_mat <- do.call(rbind, c(list(object$order), object$seasonal))
  sum_order <- colSums(order_mat)[-2]

  ar_idx <- if (sum_order[1] != 0) 1:sum_order[1] else NULL
  ma_idx <- if (sum_order[2] != 0) sum_order[1] + (1:sum_order[2]) else NULL

  sigma2 <- object$sigma2

  # Intercept, only present if the model was fitted with include.constant = TRUE
  if ("intercept" %in% names(object$coef)) {
    mu <- unname(object$coef[["intercept"]])
  } else {
    mu <- 0
  }

  msarima.sim(
    nsim = nsim,
    order = object$order,
    seasonal = object$seasonal,
    seasonal.period = object$seasonal.period,
    ar = if (length(ar_idx) > 0) object$coef[ar_idx] else NULL,
    ma = if (length(ma_idx) > 0) object$coef[ma_idx] else NULL,
    sigma2 = sigma2,
    burn = burn,
    intercept = mu
  )
}

#' Simulate from a Multi-Seasonal ARIMA Model (Direct Interface)
#'
#' @description
#' Simulates a time series from a multi-seasonal ARIMA model with specified parameters.
#' This function provides a direct interface for simulation without requiring a fitted model object.
#'
#' @param nsim Length of the simulated series (required).
#' @param order A numeric vector of length 3 specifying the non-seasonal part of the
#'   ARIMA model: the AR order (p), the degree of differencing (d), and the MA order (q).
#'   Default is \code{c(0, 0, 0)}.
#' @param seasonal A list of numeric vectors, each of length 3, specifying the seasonal
#'   components. Each vector gives the seasonal AR order (P), the seasonal differencing
#'   order (D), and the seasonal MA order (Q) of one component. Use \code{list()} for a
#'   model without seasonal components. Default is \code{list()}.
#' @param seasonal.period A numeric vector of seasonal periods, one for each element of
#'   \code{seasonal}. Default is \code{NULL}.
#' @param ar Numeric vector of AR coefficients, non-seasonal first and then one block per
#'   seasonal component. Its length must match the total AR order implied by \code{order}
#'   and \code{seasonal}. Default is \code{NULL}.
#' @param ma Numeric vector of MA coefficients, ordered as \code{ar}. Its length must match
#'   the total MA order implied by \code{order} and \code{seasonal}. Default is \code{NULL}.
#' @param sigma2 Innovation variance. Default is 1.
#' @param burn Number of burn-in observations to discard. If \code{NULL}, it is set to
#'   \code{max(100, 2 * kk)}, where \code{kk} is the largest lag in the model.
#'   Default is \code{NULL}.
#' @param intercept Constant term of the model equation, that is \eqn{\mu} in
#'   \eqn{y_t = \mu + \sum \phi_i y_{t-i} + \sum \theta_j e_{t-j} + e_t}. This is the
#'   intercept, not the mean of the process: for a stationary model the mean is
#'   \code{intercept / prod(1 - ar)}. Default is 0.
#'
#' @return A time series object of length \code{nsim}: \code{ts} with one seasonal period,
#'   \code{msts} with several.
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
  seasonal = list(),
  seasonal.period = NULL,
  ar = NULL,
  ma = NULL,
  sigma2 = 1,
  burn = NULL,
  intercept = 0
) {
  if (is.null(seasonal) || (is.list(seasonal) && length(seasonal) == 0)) {
    seasonal <- list()
    n_seasonal <- 0
    seasonal.period <- numeric(0)
  } else {
    if (!is.list(seasonal)) {
      stop("'seasonal' must be a list or NULL")
    }
    n_seasonal <- length(seasonal)

    for (i in seq_len(n_seasonal)) {
      if (!is.numeric(seasonal[[i]]) || length(seasonal[[i]]) != 3) {
        stop(sprintf("seasonal[[%d]] must be a numeric vector of length 3", i))
      }
      if (
        any(seasonal[[i]] < 0) || any(seasonal[[i]] != floor(seasonal[[i]]))
      ) {
        stop(sprintf("seasonal[[%d]] must contain non-negative integers", i))
      }
    }

    if (length(seasonal.period) != n_seasonal) {
      stop(sprintf(
        "'seasonal' and 'seasonal.period' must have the same length (%d and %d)",
        n_seasonal,
        length(seasonal.period)
      ))
    }
  }

  order_mat <- do.call(rbind, c(list(order), seasonal))
  S_vec <- c(1, seasonal.period) # Include 1 for non-seasonal

  n_ar <- sum(order_mat[, 1])
  n_ma <- sum(order_mat[, 3])

  if (n_ar > 0 && (is.null(ar) || length(ar) != n_ar)) {
    stop(sprintf(
      "'ar' must have length %d to match the AR orders in 'order' and 'seasonal'",
      n_ar
    ))
  }
  if (n_ma > 0 && (is.null(ma) || length(ma) != n_ma)) {
    stop(sprintf(
      "'ma' must have length %d to match the MA orders in 'order' and 'seasonal'",
      n_ma
    ))
  }

  sigma <- sqrt(sigma2)

  # Build AR polynomial
  if (n_ar > 0) {
    ar_comp <- rep(S_vec, order_mat[, 1])
    ar_beta <- tapply(ar, ar_comp, function(x) x, simplify = FALSE)

    poly_ar <- do.call(
      prod,
      mapply(
        function(x, y) {
          seq_par <- rep(0, length(y) * x)
          seq_par[cumsum(rep(x, length(y)))] <- y
          polynomial(c(1, -seq_par))
        },
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
        function(x, y) {
          seq_par <- rep(0, length(y) * x)
          seq_par[cumsum(rep(x, length(y)))] <- y
          polynomial(c(1, seq_par))
        },
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
  kk <- max(k_ar, k_ma, k)

  if (is.null(burn)) {
    burn <- max(100, 2 * kk)
  }

  n_total <- nsim + burn
  y <- rnorm(n_total, mean = 0, sd = sigma)
  eps <- rnorm(n_total, mean = 0, sd = sigma)

  # Start the recursion at the implied level of the process, so that the
  # burn-in does not have to absorb the whole intercept
  if (intercept != 0) {
    ar_sum <- if (length(poly) > 1) sum(as.numeric(poly)) else 1
    if (abs(ar_sum) > 1e-8) {
      y[1:kk] <- y[1:kk] + intercept / ar_sum
    }
  }

  for (j in (kk + 1):n_total) {
    y[j] <- intercept +
      sum(para * y[(j - 1):(j - k)]) +
      sum(para_ma * eps[(j - 1):(j - k_ma)]) +
      eps[j]
  }

  y_sim <- y[(burn + 1):n_total]

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
