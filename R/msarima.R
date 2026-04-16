#' Fit Multi-Seasonal ARIMA Model
#'
#' @description
#' Fits an ARIMA model with one or more seasonal components to a univariate time series.
#' This function extends the standard ARIMA framework to handle time series with
#' multiple seasonal patterns, such as daily and weekly seasonality in hourly data,
#' or weekly and yearly seasonality in daily data.
#'
#' @param y A numeric vector or time series object (ts or msts) containing the data to
#'   be modeled. If a ts object is provided and \code{seasonal.period} is not specified,
#'   the frequency from the ts object will be used as the first seasonal period. If an
#'   msts object is provided, its multiple seasonal periods will be used automatically.
#' @param order A numeric vector of length 3 specifying the non-seasonal part of the
#'   ARIMA model: the AR order (p), the degree of differencing (d), and the MA order (q).
#'   Default is \code{c(0, 0, 0)}.
#' @param seasonal A list of numeric vectors, each of length 3, specifying the seasonal
#'   components. Each vector gives the seasonal AR order (P), seasonal differencing order (D),
#'   and seasonal MA order (Q) for that seasonal component. For example,
#'   \code{list(c(1, 0, 1), c(0, 1, 1))} specifies two seasonal components.
#'   The list must contain at least one seasonal component.
#'   Default is \code{list(c(0, 0, 0), c(0, 0, 0))}.
#' @param seasonal.period A numeric vector specifying the seasonal periods corresponding
#'   to each element in \code{seasonal}. For example, \code{c(7, 30)} for weekly and
#'   monthly seasonality. If not provided and \code{y} is a ts object with frequency > 1,
#'   the frequency will be used as the first seasonal period. If \code{y} is an msts object,
#'   its seasonal periods will be used automatically. Default is \code{NULL}.
#' @param fixed A numeric vector of the same length as the total number of parameters
#'   (including sigma). If supplied, only \code{NA} entries in \code{fixed} will be varied.
#'   The order of parameters is: non-seasonal AR, seasonal AR (for each seasonal component),
#'   non-seasonal MA, seasonal MA (for each seasonal component), sigma. Use \code{fixed = NULL}
#'   to estimate all parameters, or provide a vector with \code{NA} for
#'   parameters to estimate and numeric values for parameters to fix. Default is \code{NULL}
#'   (all parameters estimated).
#' @param include.constant Should the ARIMA model include a costant term? The
#'   default is `TRUE` for undifferenced series, `FALSE` for differenced ones
#'   (where a mean would not affect the fit nor predictions).
#' @param optim_control A list of control parameters passed to \code{\link[stats]{optim}}.
#'   Default is \code{list(maxit = 200)}.
#'
#' @return An object of class \code{"msarima"}, which is a list containing:
#' \item{order}{The non-seasonal ARIMA order (p, d, q).}
#' \item{seasonal}{List of seasonal component specifications.}
#' \item{seasonal.period}{Vector of seasonal periods.}
#' \item{coef}{Named vector of estimated coefficients (AR, SAR, MA, SMA parameters).}
#' \item{var_coef}{Variance-covariance matrix of the coefficient estimates.}
#' \item{mask}{Logical vector indicating which parameters were fixed.}
#' \item{loglik}{Log-likelihood of the fitted model.}
#' \item{aic}{Akaike Information Criterion.}
#' \item{bic}{Bayesian Information Criterion.}
#' \item{aicc}{Corrected Akaike Information Criterion.}
#' \item{sigma2}{Estimated variance of the residuals.}
#' \item{x}{The original time series.}
#' \item{fitted}{Fitted values.}
#' \item{residuals}{Model residuals.}
#' \item{method}{Description of the fitted model.}
#' \item{series}{Name of the series.}
#' \item{nobs}{Number of observations in the original series.}
#' \item{n.used}{Number of observations used in estimation (after differencing).}
#' \item{call}{The matched call.}
#' \item{convergence}{Convergence code from \code{optim} (0 indicates successful convergence).}
#'
#' @details
#' The model is estimated using maximum likelihood via numerical optimization with
#' \code{\link[stats]{optim}} using the BFGS algorithm.
#'
#' Standard errors are computed from the inverse of the observed information matrix
#' (negative Hessian). If the Hessian is singular, a warning is issued and standard
#' errors are set to \code{NA}.
#'
#' @note
#' \itemize{
#'   \item The series must have at least 10 observations.
#'   \item The \code{seasonal} list must contain at least one seasonal component (it cannot be empty).
#'   \item Seasonal periods must be at least 2 and less than half the series length.
#'   \item Differencing orders (d, D for each seasonal component) are restricted to 0 or 1.
#'   \item The function may take longer to converge for complex models with many parameters.
#'   \item If input is ts or msts, the output x, fitted, and residuals preserve the same structure.
#'   \item If input is numeric and seasonal periods are specified, output is msts.
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
#' # Example 2: Using msts object (multiple periods detected automatically)
#' y_msts <- msts(y, seasonal.periods = c(4, 11))
#' fit_msts <- msarima(y_msts, order = c(1, 0, 0),
#'                     seasonal = list(c(1, 0, 0), c(0, 0, 1)))
#' # Both periods (4 and 11) are used automatically
#'
#' # Example 4: Fixing parameters
#' # Note: fixed vector length must match total parameters (AR + SAR + MA + SMA + sigma)
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
  seasonal = list(c(0, 0, 0), c(0, 0, 0)),
  seasonal.period = NULL,
  fixed = NULL,
  include.constant = FALSE,
  optim_control = list(maxit = 200)
) {
  # Input validation
  if (!is.numeric(y)) {
    stop("'y' must be numeric")
  }
  if (any(is.na(y))) {
    stop("'y' contains missing values")
  }
  if (length(y) < 10) {
    stop("Series too short")
  }

  series <- deparse1(substitute(y))

  # Extract frequency from ts or msts object if available
  y_frequency <- NULL
  if (inherits(y, "msts")) {
    y_frequency <- attr(y, "msts")
  } else if (is.ts(y)) {
    y_frequency <- frequency(y)
  }

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

  # Handle seasonal.period
  if (is.null(seasonal.period)) {
    seasonal.period <- y_frequency
  }

  if (length(seasonal.period) != n_seasonal) {
    stop(
      "'seasonal.period' must have the same lenght of 'seasonal.period'"
    )
  }

  # Store seasonal orders for easy access
  seasonal_orders <- seasonal

  n <- length(y)

  # Combine orders for easier handling (P, D, Q)
  order_mat <- do.call(rbind, c(list(order), seasonal_orders))
  S_vec <- seasonal.period

  if (any(order_mat[, 2] > 0) && include.constant) {
    warning("No constant term fitted as the order of difference is 1 or more.")
    include.constant <- FALSE
  }

  # Validate seasonal periods
  for (s in 1:n_seasonal) {
    if (!is.na(S_vec[s]) && (S_vec[s] < 2 || S_vec[s] > n / 2)) {
      stop(sprintf("Invalid seasonal period S%d", s))
    }
  }

  if (any(order_mat[, 2] > 1)) {
    stop("Differencing orders must be 0 or 1")
  }

  # Apply differencing
  xx <- as.numeric(y)
  nstar <- length(xx)

  # Total number of parameters (including sigma)
  n_params <- sum(order_mat[, c(1, 3)]) # AR and MA parameters only
  n_params_total <- n_params + 1 + include.constant # +1 for sigma, +1 mu

  # Handle fixed parameters
  if (!is.null(fixed)) {
    if (length(fixed) != n_params_total) {
      stop(sprintf(
        "'fixed' must be NULL or a vector of length %d (total number of parameters including sigma)",
        n_params_total
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

  # Log-likelihood function (vectorized version)
  loglik <- function(beta_free) {
    # Reconstruct full parameter vector (including sigma)
    beta <- numeric(n_params_total)

    if (all_fixed) {
      beta <- fixed
    } else {
      beta[mask] <- fixed[mask]
      beta[!mask] <- beta_free
    }

    # Extract sigma (last element)
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
          function(x, y) polynomial(c(1, rep(0, x - 1), -y)),
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
          function(x, y) polynomial(c(1, rep(0, x - 1), y)),
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
          function(x, y) polynomial(c(1, rep(0, x - 1), y)),
          x = c(1, seasonal.period),
          y = order_mat[, 2],
          SIMPLIFY = FALSE
        )
      )
    } else {
      poly_i <- polynomial(1)
    }

    # Construct total polynomials
    poly <- poly_ar * poly_i

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
    #xx_start <- rep(mean(xx[1:20]), kk)
    xxl <- c(xx_start, xx)
    nxxl <- length(xxl)
    epsl <- rep(0, nxxl)

    # Residual calculation
    for (j in (kk + 1):nxxl) {
      epsl[j] <- xxl[j] -
        mu -
        (sum(para_ar * xxl[(j - 1):(j - k_ar)]) +
          sum(para_ma * epsl[(j - 1):(j - k_ma)]))
    }

    # Remove burn-in
    residuals <- epsl[(kk + 1):nxxl]

    # Vectorized log-likelihood calculation
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
  beta_ini <- beta_ini * 0.1
  beta_ini <- as.vector(beta_ini)
  beta_ini <- beta_ini[!is.na(beta_ini)]
  beta_ini_full <- c(beta_ini, sd(xx))

  if (include.constant) {
    beta_ini_full <- c(beta_ini_full, mean(xx))
  }

  # Use fixed values where specified
  if (any(mask)) {
    beta_ini_full[mask] <- fixed[mask]
  }

  # Extract only free parameters for optimization
  beta_ini <- beta_ini_full[!mask]

  # Optimize log-likelihood
  if (all_fixed) {
    fit <- list(
      par = fixed,
      value = loglik(numeric(0)),
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

  # Check convergence
  if (fit$convergence != 0) {
    warning(paste("Optimization did not converge. Code:", fit$convergence))
  }

  # Reconstruct full parameter vector
  if (all_fixed) {
    beta_full <- fixed
  } else {
    beta_full <- numeric(n_params_total)
    beta_full[mask] <- fixed[mask]
    beta_full[!mask] <- fit$par
  }

  LLO <- loglik(beta_full)

  # Extract estimated parameters
  sigma2_hat <- beta_full[n_params_total - include.constant]^2
  if (include.constant) {
    mu_hat <- beta_full[n_params_total]
  }

  # Coefficients (AR and MA parameters only, exclude sigma)
  Coef <- beta_full[include.constant - n_params_total]

  # Calculate variance-covariance matrix
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
      warning("Hessian is singular. Standard errors may be unreliable.")
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
      include.constant - n_params_total
    ]

    if (any(is.na(fit$hessian)) || det(-fit$hessian) <= 0) {
      warning("Hessian is singular. Standard errors may be unreliable.")
      var_coef <- matrix(NA, nrow = n_params, ncol = n_params)
    } else {
      j0 <- (-fit$hessian) / nstar
      var_coef <- solve(j0) / nstar
    }
  }

  # Calculate information criteria
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
        "AR"
      } else {
        paste0("SAR", i - 1, ".", 1:order_mat[i, 1])
      }
    } else {
      NA
    }
  })
  ar_name <- ar_name[!is.na(ar_name)]

  ma_name <- sapply(1:NROW(order_mat), function(i) {
    if (order_mat[i, 3] > 0) {
      if (i == 1) {
        "MA"
      } else {
        paste0("SMA", i - 1, ".", 1:order_mat[i, 3])
      }
    } else {
      NA
    }
  })
  ma_name <- ma_name[!is.na(ma_name)]
  coef_names <- c(ar_name, ma_name)

  if (include.constant) {
    coef_names <- c(coef_names, "intercept")
  }

  resid <- LLO$residuals
  fits <- xx - LLO$residuals

  if (is.null(tsp(y))) {
    valid_periods <- seasonal.period[!is.na(seasonal.period)]
    if (length(valid_periods) > 0) {
      y <- msts(y, seasonal.periods = valid_periods)
    } else {
      y <- ts(y)
    }
  }

  # Create method description
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

  # Create output object
  out <- structure(
    list(
      # Model specification
      order = order,
      seasonal = seasonal,
      seasonal.period = seasonal.period,

      # Coefficients
      coef = setNames(Coef, coef_names),
      var_coef = var_coef,
      mask = mask_params,

      # Model fit statistics
      loglik = loglik_val,
      aic = aic,
      bic = bic,
      aicc = aicc,
      sigma2 = sigma2_hat,

      # Series and residuals
      x = y,
      fitted = copy_msts(y, fits),
      residuals = copy_msts(y, resid),

      # Additional info
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
