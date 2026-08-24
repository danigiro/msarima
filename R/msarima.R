#' Fit a Multi-Seasonal ARIMA Model
#'
#' @description
#' Fits an ARIMA model with one or more seasonal components to a univariate time series.
#' This function extends the standard ARIMA framework to time series with several
#' seasonal patterns, such as daily and weekly seasonality in hourly data, or weekly
#' and yearly seasonality in daily data.
#'
#' @param y A numeric vector or time series object (\code{ts} or \code{msts}) containing
#'   the data to be modelled. If \code{y} is a \code{ts} object and \code{seasonal.period}
#'   is not supplied, its frequency is used as the first seasonal period; if \code{y} is
#'   an \code{msts} object, its seasonal periods are used.
#' @param order A numeric vector of length 3 specifying the non-seasonal part of the
#'   ARIMA model: the AR order (p), the degree of differencing (d), and the MA order (q).
#'   Default is \code{c(0, 0, 0)}.
#' @param seasonal A list of numeric vectors, each of length 3, specifying the seasonal
#'   components. Each vector gives the seasonal AR order (P), the seasonal differencing
#'   order (D), and the seasonal MA order (Q) of one component. For example,
#'   \code{list(c(1, 0, 1), c(0, 1, 1))} specifies two seasonal components. Use
#'   \code{list()} for a model without seasonal components. Default is \code{list()}.
#' @param seasonal.period A numeric vector of seasonal periods, one for each element of
#'   \code{seasonal}, for example \code{c(7, 30)} for weekly and monthly seasonality.
#'   If \code{NULL}, the periods are taken from \code{y} when it is a \code{ts} or
#'   \code{msts} object. Default is \code{NULL}.
#' @param fixed A numeric vector with one entry per parameter, used to hold some of them
#'   at a given value: only the \code{NA} entries are estimated. The parameters are
#'   ordered as non-seasonal AR, seasonal AR (one block per seasonal component),
#'   non-seasonal MA, seasonal MA (one block per seasonal component), sigma, and the
#'   constant term when \code{include.constant = TRUE}. Default is \code{NULL}, meaning
#'   that all parameters are estimated.
#' @param include.constant Should the model include a constant term? The constant is the
#'   intercept of the model equation, not the mean of the series: for a stationary model
#'   the two are related by \code{mean = intercept / prod(1 - ar)}. A constant cannot be
#'   combined with differencing and is dropped, with a warning, when \code{d} or any
#'   \code{D} is non-zero. Default is \code{FALSE}.
#' @param optim_control A list of control parameters passed to \code{\link[stats]{optim}}.
#'   Default is \code{list(maxit = 200)}.
#'
#' @return An object of class \code{"msarima"}, which is a list containing:
#' \item{order}{The non-seasonal ARIMA order (p, d, q).}
#' \item{seasonal}{List of seasonal component specifications.}
#' \item{seasonal.period}{Vector of seasonal periods.}
#' \item{coef}{Named vector of estimated coefficients: the AR, SAR, MA and SMA
#'   parameters, followed by the intercept if \code{include.constant = TRUE}.}
#' \item{var_coef}{Variance-covariance matrix of the coefficient estimates.}
#' \item{mask}{Logical vector indicating which parameters were fixed.}
#' \item{loglik}{Log-likelihood of the fitted model.}
#' \item{aic}{Akaike Information Criterion.}
#' \item{bic}{Bayesian Information Criterion.}
#' \item{aicc}{Corrected Akaike Information Criterion.}
#' \item{sigma2}{Estimated innovation variance.}
#' \item{x}{The original time series.}
#' \item{fitted}{One-step-ahead fitted values, on the scale of the original series.}
#' \item{residuals}{Model residuals. Together with \code{fitted}, the first
#'   \code{d + sum(D * S)} values are \code{NA}, since differencing leaves no
#'   observation to predict there.}
#' \item{method}{Description of the fitted model.}
#' \item{series}{Name of the series.}
#' \item{nobs}{Number of observations in the original series.}
#' \item{n.used}{Number of observations used in estimation (after differencing).}
#' \item{call}{The matched call.}
#' \item{convergence}{Convergence code from \code{optim} (0 indicates successful convergence).}
#'
#' @details
#' The model is estimated by conditional maximum likelihood: the residuals are computed
#' recursively and the resulting Gaussian log-likelihood is maximised numerically with
#' \code{\link[stats]{optim}}, using the BFGS algorithm.
#'
#' When \code{d} or any \code{D} is non-zero, the series is differenced before estimation
#' and the ARMA part is estimated on the differenced data. The first
#' \code{d + sum(D * S)} observations are therefore not used, which is why the residuals
#' start with as many \code{NA} values.
#'
#' Standard errors are computed from the inverse of the observed information matrix
#' (negative Hessian). If the Hessian is singular, a warning is issued and the standard
#' errors are set to \code{NA}.
#'
#' @note
#' \itemize{
#'   \item The series must have at least 10 observations, both before and after differencing.
#'   \item Seasonal periods must be at least 2 and at most half the length of the series.
#'   \item If \code{y} is a \code{ts} or \code{msts} object, \code{x}, \code{fitted} and
#'     \code{residuals} keep the same time series structure. If \code{y} is a plain numeric
#'     vector, they are returned as \code{ts} with one seasonal period, or as \code{msts}
#'     with several.
#' }
#'
#' @examples
#' # Example 1: ARIMA with two seasonal components
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
#' print(fit)
#' summary(fit)
#'
#' # Example 2: msts input, seasonal periods detected automatically
#' y_msts <- msts(y, seasonal.periods = c(4, 11))
#' fit_msts <- msarima(y_msts, order = c(1, 0, 0),
#'                     seasonal = list(c(1, 0, 0), c(0, 0, 1)))
#'
#' # Example 3: differencing
#' fit_diff <- msarima(cumsum(y), order = c(1, 1, 0),
#'                     seasonal = list(c(0, 0, 1)),
#'                     seasonal.period = 4)
#' head(residuals(fit_diff))  # the first residual is NA
#'
#' # Example 4: fixing parameters
#' # Note: the length of 'fixed' must match AR + SAR + MA + SMA + sigma
#' fit_fixed <- msarima(y, order = c(1, 0, 0),
#'                      seasonal = list(c(1, 0, 0), c(0, 0, 1)),
#'                      seasonal.period = c(4, 11),
#'                      fixed = c(0.5, NA, NA, NA))
#' AIC(fit); AIC(fit_fixed)
#'
#' @export
msarima <- function(
  y,
  order = c(0, 0, 0),
  seasonal = list(),
  seasonal.period = NULL,
  fixed = NULL,
  include.constant = FALSE,
  optim_control = list(maxit = 200)
) {
  if (!is.numeric(y)) {
    stop("'y' must be numeric")
  }
  if (any(is.na(y))) {
    stop("'y' contains missing values")
  }
  if (length(y) < 10) {
    stop(sprintf("'y' has %d observations, at least 10 are needed", length(y)))
  }

  series <- deparse1(substitute(y))

  y_frequency <- NULL
  if (inherits(y, "msts")) {
    y_frequency <- attr(y, "msts")
  } else if (is.ts(y)) {
    y_frequency <- frequency(y)
  }

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

    if (is.null(seasonal.period)) {
      seasonal.period <- y_frequency
    }

    if (length(seasonal.period) != n_seasonal) {
      stop(sprintf(
        "'seasonal' and 'seasonal.period' must have the same length (%d and %d)",
        n_seasonal,
        length(seasonal.period)
      ))
    }
  }

  seasonal_orders <- seasonal

  n <- length(y)

  order_mat <- do.call(rbind, c(list(order), seasonal_orders))
  S_vec <- seasonal.period

  if (any(order_mat[, 2] > 0) && include.constant) {
    warning("No constant term fitted: the model includes differencing")
    include.constant <- FALSE
  }

  if (n_seasonal > 0) {
    for (s in 1:n_seasonal) {
      if (!is.na(S_vec[s]) && (S_vec[s] < 2 || S_vec[s] > n / 2)) {
        stop(sprintf(
          "Invalid 'seasonal.period[%d]' = %g: seasonal periods must be between 2 and %g (half the length of 'y')",
          s,
          S_vec[s],
          n / 2
        ))
      }
    }
  }

  # Apply differencing
  y_num <- as.numeric(y)
  xx <- y_num
  S_all <- c(1, S_vec)
  for (i in seq_len(NROW(order_mat))) {
    if (order_mat[i, 2] > 0) {
      xx <- diff(xx, lag = S_all[i], differences = order_mat[i, 2])
    }
  }
  n_diff <- n - length(xx) # observations lost to differencing
  nstar <- length(xx)
  if (nstar < 10) {
    stop(sprintf(
      "Only %d observations are left after differencing, at least 10 are needed",
      nstar
    ))
  }

  n_params <- sum(order_mat[, c(1, 3)]) # AR and MA parameters only
  n_params_total <- n_params + 1 + include.constant # +1 for sigma, +1 mu

  if (!is.null(fixed)) {
    if (length(fixed) != n_params_total) {
      stop(sprintf(
        "'fixed' has length %d, but the model has %d parameters: %d AR/MA coefficients, sigma%s",
        length(fixed),
        n_params_total,
        n_params,
        if (include.constant) " and the constant" else ""
      ))
    }
    mask <- !is.na(fixed)
    sigma_fixed <- mask[n_params_total]
    mask_params <- mask[1:n_params]
  } else {
    fixed <- rep(NA, n_params_total)
    mask <- rep(FALSE, n_params_total)
    sigma_fixed <- FALSE
    mask_params <- rep(FALSE, n_params)
  }

  all_fixed <- all(mask)

  loglik <- function(beta_free) {
    beta <- numeric(n_params_total)

    if (all_fixed) {
      beta <- fixed
    } else {
      beta[mask] <- fixed[mask]
      beta[!mask] <- beta_free
    }

    sigma <- beta[n_params_total - include.constant]

    if (include.constant) {
      mu <- beta[n_params_total]
    } else {
      mu <- 0
    }

    if (sum(order_mat[, 1]) != 0) {
      ar_comp <- rep(c(1, seasonal.period), order_mat[, 1])
      ar_beta <- beta[1:length(ar_comp)]
      ar_beta <- tapply(ar_beta, ar_comp, function(x) x, simplify = FALSE)

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
      ar_comp <- NULL
      poly_ar <- polynomial(1)
    }

    if (sum(order_mat[, 3]) != 0) {
      ma_comp <- rep(c(1, seasonal.period), order_mat[, 3])
      ma_beta <- beta[length(ar_comp) + 1:length(ma_comp)]
      ma_beta <- tapply(ma_beta, ma_comp, function(x) x, simplify = FALSE)
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

    if (sum(order_mat[, 2]) != 0) {
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
          x = c(1, seasonal.period),
          y = order_mat[, 2],
          SIMPLIFY = FALSE
        )
      )
    } else {
      poly_i <- polynomial(1)
    }

    # Construct total polynomials. `xx` has already been differenced, so the
    # recursion only uses the stationary AR polynomial (poly_i is kept for
    # reporting / forecasting).
    poly <- poly_ar

    if (length(poly) > 1) {
      para_ar <- -poly[2:length(poly)]
    } else {
      para_ar <- 0
    }

    if (length(poly_ma) > 1) {
      para_ma <- poly_ma[2:length(poly_ma)]
    } else {
      para_ma <- 0
    }

    k_ar <- length(para_ar)
    k_ma <- length(para_ma)
    nxx <- length(xx)
    kk <- max(k_ar, k_ma)
    xx_start <- rep(mean(xx[1:min(20, nxx / 5)]), kk)
    xxl <- c(xx_start, xx)
    nxxl <- length(xxl)
    epsl <- rep(0, nxxl)

    for (j in (kk + 1):nxxl) {
      epsl[j] <- xxl[j] -
        mu -
        (sum(para_ar * xxl[(j - 1):(j - k_ar)]) +
          sum(para_ma * epsl[(j - 1):(j - k_ma)]))
    }

    residuals <- epsl[(kk + 1):nxxl]

    sigma2 <- sigma^2
    n_eff <- nxxl - kk
    ss <- sum(residuals^2)
    l <- -0.5 * n_eff * log(sigma2) - 0.5 * ss / sigma2
    return(list(
      value = l,
      residuals = residuals,
      poly_ar = poly_ar,
      poly_i = poly_i,
      poly_ma = poly_ma
    ))
  }

  beta_ini <- order_mat[, -2]
  beta_ini[beta_ini == 0] <- NA
  beta_ini <- rep(0.1, sum(beta_ini, na.rm = TRUE))
  beta_ini <- beta_ini[!is.na(beta_ini)]
  beta_ini_full <- c(beta_ini, sd(xx))

  if (include.constant) {
    beta_ini_full <- c(beta_ini_full, mean(xx))
  }

  if (any(mask)) {
    beta_ini_full[mask] <- fixed[mask]
  }

  beta_ini <- beta_ini_full[!mask]

  if (all_fixed) {
    fit <- list(
      par = fixed,
      value = loglik(numeric(0))$value,
      convergence = 0,
      hessian = matrix(0, nrow = 0, ncol = 0)
    )
  } else {
    fit <- optim(
      par = beta_ini,
      fn = function(b) loglik(b)$value,
      gr = NULL,
      method = "BFGS",
      control = c(list(fnscale = -1), optim_control),
      hessian = TRUE
    )
  }

  if (fit$convergence != 0) {
    if (fit$convergence == 1) {
      warning(
        "Optimization did not converge: iteration limit reached, try increasing 'maxit' in 'optim_control'"
      )
    } else {
      warning(sprintf(
        "Optimization did not converge (optim code %d), the estimates may not be reliable",
        fit$convergence
      ))
    }
  }

  if (all_fixed) {
    beta_full <- fixed
  } else {
    beta_full <- numeric(n_params_total)
    beta_full[mask] <- fixed[mask]
    beta_full[!mask] <- fit$par
  }

  # loglik() expects only the free parameters, not the full vector
  LLO <- loglik(beta_full[!mask])

  sigma2_hat <- beta_full[n_params_total - include.constant]^2
  if (include.constant) {
    mu_hat <- beta_full[n_params_total]
  }

  if (all_fixed) {
    var_coef <- matrix(0, nrow = n_params, ncol = n_params)
  } else if (any(mask_params)) {
    hess_full <- matrix(0, nrow = n_params_total, ncol = n_params_total)
    free_idx <- which(!mask)
    if (length(free_idx) > 0) {
      hess_free <- fit$hessian[1:length(free_idx), 1:length(free_idx)]
      hess_full[free_idx, free_idx] <- hess_free
    }
    fit$hessian <- hess_full[1:n_params, 1:n_params]

    free_params_idx <- which(!mask_params)
    if (
      length(free_params_idx) > 0 &&
        (any(is.na(fit$hessian[free_params_idx, free_params_idx])) ||
          det(-fit$hessian[free_params_idx, free_params_idx]) <= 0)
    ) {
      warning("Hessian is singular: standard errors are not available")
      var_coef <- matrix(NA, nrow = n_params, ncol = n_params)
    } else if (length(free_params_idx) > 0) {
      j0_free <- (-fit$hessian[free_params_idx, free_params_idx]) / nstar
      var.coef_free <- solve(j0_free) / nstar
      var_coef <- matrix(NA, nrow = n_params, ncol = n_params)
      var_coef[free_params_idx, free_params_idx] <- var.coef_free
    } else {
      var_coef <- matrix(NA, nrow = n_params, ncol = n_params)
    }
  } else {
    fit$hessian <- fit$hessian[
      include.constant - n_params_total,
      include.constant - n_params_total,
      drop = FALSE
    ]
    if (any(dim(fit$hessian) == 0)) {
      var_coef <- NULL
    } else if (any(is.na(fit$hessian)) || det(-fit$hessian) <= 0) {
      warning("Hessian is singular: standard errors are not available")
      var_coef <- matrix(NA, nrow = n_params, ncol = n_params)
    } else {
      j0 <- (-fit$hessian) / nstar
      var_coef <- solve(j0) / nstar
    }
  }

  if (all_fixed) {
    npar_total <- 0
  } else {
    npar_total <- sum(!mask)
  }
  loglik_val <- fit$value
  aic <- -2 * loglik_val + 2 * npar_total
  bic <- -2 * loglik_val + log(nstar) * npar_total
  if (npar_total > 0) {
    aicc <- aic + (2 * npar_total * (npar_total + 1)) / (nstar - npar_total - 1)
  } else {
    aicc <- aic
  }

  ar_name <- sapply(1:NROW(order_mat), function(i) {
    if (order_mat[i, 1] > 0) {
      if (i == 1) {
        paste0("AR(", 1:order_mat[i, 1], ")")
      } else {
        paste0("SAR(", i - 1, ",", 1:order_mat[i, 1], ")")
      }
    } else {
      NA
    }
  })
  ar_name <- unlist(ar_name[!is.na(ar_name)])

  ma_name <- sapply(1:NROW(order_mat), function(i) {
    if (order_mat[i, 3] > 0) {
      if (i == 1) {
        paste0("MA(", 1:order_mat[i, 3], ")")
      } else {
        paste0("SMA(", i - 1, ",", 1:order_mat[i, 3], ")")
      }
    } else {
      NA
    }
  })
  ma_name <- unlist(ma_name[!is.na(ma_name)])
  coef_names <- c(ar_name, ma_name)

  if (include.constant) {
    coef_names <- c(coef_names, "intercept")
  }

  resid <- LLO$residuals
  fits <- y_num[(n_diff + 1):n] - resid

  if (is.null(tsp(y))) {
    valid_periods <- seasonal.period[!is.na(seasonal.period)]
    if (length(valid_periods) > 0) {
      y <- msts(y, seasonal.periods = valid_periods)
    } else {
      y <- ts(y)
    }
  }

  method_desc <- sprintf(
    "MSARIMA(%d,%d,%d)",
    order_mat[1, 1],
    order_mat[1, 2],
    order_mat[1, 3]
  )

  method_desc <- paste0(
    c(
      method_desc,
      sprintf(
        "(%d,%d,%d)[%d]",
        order_mat[-1, 1],
        order_mat[-1, 2],
        order_mat[-1, 3],
        seasonal.period
      )
    ),
    collapse = ""
  )

  coef_vec <- beta_full[include.constant - n_params_total]
  if (length(coef_vec) == 0) {
    coef_vec <- NULL
  } else {
    coef_vec <- setNames(coef_vec, coef_names)
  }

  out <- structure(
    list(
      order = order,
      seasonal = seasonal,
      seasonal.period = seasonal.period,

      coef = coef_vec,
      var_coef = var_coef,
      mask = mask_params,

      loglik = loglik_val,
      aic = aic,
      bic = bic,
      aicc = aicc,
      sigma2 = sigma2_hat,

      x = y,
      fitted = copy_msts(y, fits),
      residuals = copy_msts(y, resid),

      method = method_desc,
      series = series,
      nobs = n,
      n.used = nstar,
      call = match.call(),
      convergence = fit$convergence
    ),
    class = "msarima"
  )
  return(out)
}
