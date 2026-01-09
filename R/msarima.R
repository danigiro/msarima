#' Fit Multi-Seasonal ARIMA Model
#'
#' @description
#' Fits an ARIMA model with up to two seasonal components to a univariate time series.
#' This function extends the standard ARIMA framework to handle time series with
#' multiple seasonal patterns, such as daily and weekly seasonality in hourly data,
#' or weekly and yearly seasonality in daily data.
#'
#' @param y A numeric vector or time series object containing the data to be modeled.
#' @param order A numeric vector of length 3 specifying the non-seasonal part of the
#'   ARIMA model: the AR order (p), the degree of differencing (d), and the MA order (q).
#'   Default is \code{c(0, 0, 0)}.
#' @param seasonal A list with components \code{order} and \code{period} specifying
#'   the first seasonal part of the model: \code{order} is a numeric vector of length 3
#'   giving the seasonal AR order (P), seasonal differencing order (D), and seasonal MA
#'   order (Q); \code{period} is the seasonal period.
#' @param seasonal2 A list with components \code{order} and \code{period} specifying
#'   the second seasonal part of the model: same structure as \code{seasonal}.
#' @param optim_control A list of control parameters passed to \code{\link[stats]{optim}}.
#'   Default is \code{list(maxit = 200)}.
#'
#' @return An object of class \code{"msarima"}, which is a list containing:
#' \item{order}{The non-seasonal ARIMA order (p, d, q).}
#' \item{seasonal}{List containing the first seasonal component specification.}
#' \item{seasonal2}{List containing the second seasonal component specification.}
#' \item{coef}{Named vector of estimated coefficients (AR, SAR, MA, SMA parameters).}
#' \item{var.coef}{Variance-covariance matrix of the coefficient estimates.}
#' \item{loglik}{Log-likelihood of the fitted model.}
#' \item{aic}{Akaike Information Criterion.}
#' \item{bic}{Bayesian Information Criterion.}
#' \item{aicc}{Corrected Akaike Information Criterion.}
#' \item{sigma2}{Estimated variance of the residuals.}
#' \item{x}{The original time series.}
#' \item{fitted}{Fitted values.}
#' \item{residuals}{Model residuals.}
#' \item{series}{Name of the series.}
#' \item{method}{Estimation method used.}
#' \item{nobs}{Number of observations in the original series.}
#' \item{n.used}{Number of observations used in estimation (after differencing).}
#' \item{call}{The matched call.}
#' \item{convergence}{Convergence code from \code{optim} (0 indicates successful convergence).}
#'
#' @details
#' The function fits a multi-seasonal ARIMA model of the form:
#'
#' \deqn{\Phi_p(B)\Phi_{P1}(B^{s1})\Phi_{P2}(B^{s2})(1-B)^d(1-B^{s1})^{D1}(1-B^{s2})^{D2}y_t = \Theta_q(B)\Theta_{Q1}(B^{s1})\Theta_{Q2}(B^{s2})\epsilon_t}
#'
#' where:
#' \itemize{
#'   \item \eqn{\Phi_p(B)} is the non-seasonal AR polynomial of order p
#'   \item \eqn{\Phi_{P1}(B^{s1})} is the first seasonal AR polynomial of order P1 with period s1
#'   \item \eqn{\Phi_{P2}(B^{s2})} is the second seasonal AR polynomial of order P2 with period s2
#'   \item \eqn{\Theta_q(B)} is the non-seasonal MA polynomial of order q
#'   \item \eqn{\Theta_{Q1}(B^{s1})} is the first seasonal MA polynomial of order Q1
#'   \item \eqn{\Theta_{Q2}(B^{s2})} is the second seasonal MA polynomial of order Q2
#'   \item d, D1, D2 are the differencing orders (must be 0 or 1)
#'   \item \eqn{\epsilon_t} is white noise with variance \eqn{\sigma^2}
#' }
#'
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
#'   \item Seasonal periods must be at least 2 and less than half the series length.
#'   \item Differencing orders (d, D1, D2) are restricted to 0 or 1.
#'   \item The function may take longer to converge for complex models with many parameters.
#'   \item For very long series or high-order models, consider using smaller \code{optim.control$maxit}.
#' }
#'
#' @examples
#' # Simulate data with two seasonal patterns
#' set.seed(123)
#' n <- 500
#' t <- 1:n
#' y <- 10 + 0.05*t + 3*sin(2*pi*t/7) + 2*sin(2*pi*t/30) + rnorm(n, sd=0.5)
#'
#' # Fit ARIMA(1,0,1) with weekly seasonality (period=7)
#' fit1 <- msarima(y,
#'                 order = c(1, 0, 1),
#'                 seasonal = list(order = c(1, 0, 1), period = 7))
#' print(fit1)
#' summary(fit1)
#'
#' # Fit model with two seasonal components
#' fit2 <- msarima(y,
#'                 order = c(1, 1, 1),
#'                 seasonal = list(order = c(1, 0, 1), period = 7),
#'                 seasonal2 = list(order = c(0, 1, 1), period = 30))
#' print(fit2)
#'
#' # Extract components
#' coef(fit2)
#' fitted(fit2)
#' residuals(fit2)
#' AIC(fit2)
#' BIC(fit2)
#'
#' # Plot residuals
#' plot(residuals(fit2), type = "l", main = "Residuals")
#' acf(residuals(fit2), main = "ACF of Residuals")
#'
#' @export
msarima <- function(
  y,
  order = c(0, 0, 0),
  seasonal = list(order = c(0, 0, 0), period = NA),
  seasonal2 = list(order = c(0, 0, 0), period = NA),
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

  check_seasonal(seasonal, "seasonal")
  check_seasonal(seasonal2, "seasonal2")

  # Store original series
  origdata <- y
  n <- length(y)

  # Parse orders
  p <- order[1]
  di <- order[2]
  q <- order[3]
  P1 <- seasonal$order[1]
  Di1 <- seasonal$order[2]
  Q1 <- seasonal$order[3]
  P2 <- seasonal2$order[1]
  Di2 <- seasonal2$order[2]
  Q2 <- seasonal2$order[3]
  S1 <- seasonal$period
  S2 <- seasonal2$period

  # Validate seasonal periods
  if (!is.na(S1) && (S1 < 2 || S1 > n / 2)) {
    stop("Invalid seasonal period S1")
  }
  if (!is.na(S2) && (S2 < 2 || S2 > n / 2)) {
    stop("Invalid seasonal period S2")
  }

  # Apply differencing
  xx <- origdata
  if (di == 1) {
    xx <- diff(xx, lag = 1)
  }
  if (Di1 == 1) {
    xx <- diff(xx, lag = S1)
  }
  if (Di2 == 1) {
    xx <- diff(xx, lag = S2)
  }

  nstar <- length(xx)

  # Log-likelihood function for the stationary series
  loglik <- function(beta) {
    n_beta <- length(beta)
    sigma <- beta[n_beta]

    # Build AR polynomial components
    poly.ar <- if (p == 0) polynomial(1) else polynomial(c(1, -beta[1:p]))

    poly.ar.s1 <- polynomial(1)
    if (P1 > 0) {
      for (j in 1:P1) {
        poly.ar.s1 <- polynomial(c(poly.ar.s1, rep(0, S1 - 1), -beta[p + j]))
      }
    }

    poly.ar.s2 <- polynomial(1)
    if (P2 > 0) {
      for (j in 1:P2) {
        poly.ar.s2 <- polynomial(c(
          poly.ar.s2,
          rep(0, S2 - 1),
          -beta[p + P1 + j]
        ))
      }
    }

    # Build MA polynomial components
    poly.ma <- if (q == 0) {
      polynomial(1)
    } else {
      polynomial(c(1, beta[(p + P1 + P2 + 1):(p + P1 + P2 + q)]))
    }

    poly.ma.s1 <- polynomial(1)
    if (Q1 > 0) {
      for (j in 1:Q1) {
        poly.ma.s1 <- polynomial(c(
          poly.ma.s1,
          rep(0, S1 - 1),
          beta[(p + P1 + P2 + q + j)]
        ))
      }
    }

    poly.ma.s2 <- polynomial(1)
    if (Q2 > 0) {
      for (j in 1:Q2) {
        poly.ma.s2 <- polynomial(c(
          poly.ma.s2,
          rep(0, S2 - 1),
          beta[(p + P1 + P2 + q + Q1 + j)]
        ))
      }
    }

    # Build integration polynomials (assumes d, D1, D2 can only be 0 or 1)
    poly.i <- polynomial(c(1, -di))
    if (!is.na(S1)) {
      poly.s1.i <- polynomial(c(1, rep(0, S1 - 1), -Di1))
    } else {
      poly.s1.i <- polynomial(1)
    }

    if (!is.na(S2)) {
      poly.s2.i <- polynomial(c(1, rep(0, S2 - 1), -Di2))
    } else {
      poly.s2.i <- polynomial(1)
    }

    # Construct total AR and MA polynomials
    poly.ar.tot <- poly.ar * poly.ar.s1 * poly.ar.s2
    poly.ma.tot <- poly.ma * poly.ma.s1 * poly.ma.s2
    poly <- poly.ar * poly.ar.s1 * poly.ar.s2 * poly.i * poly.s1.i * poly.s2.i

    # Extract parameters (set to 0 if all orders are 0)
    para.ar <- if ((p + P1 + P2) == 0) {
      0
    } else {
      -poly.ar.tot[2:length(poly.ar.tot)]
    }
    para.ma <- if ((q + Q1 + Q2) == 0) 0 else poly.ma.tot[2:length(poly.ma.tot)]
    para <- if ((p + P1 + P2 + di + Di1 + Di2) == 0) {
      0
    } else {
      -poly[2:length(poly)]
    }

    k.ar <- length(para.ar)
    k.ma <- length(para.ma)
    k <- length(para)

    # Initialize extended series for residual calculation
    nxx <- length(xx)
    kk <- max(k.ar, k.ma)
    xx.start <- rep(mean(xx[1:20]), kk)
    xxl <- c(xx.start, xx)
    nxxl <- length(xxl)
    epsl <- rep(0, nxxl)

    # Calculate residuals (epsilon series)
    for (j in (kk + 1):nxxl) {
      epsl[j] <- xxl[j] -
        (sum(para.ar * xxl[(j - 1):(j - k.ar)]) +
          sum(para.ma * epsl[(j - 1):(j - k.ma)]))
    }

    # Calculate conditional means
    mu <- rep(0, nxxl)
    for (j in (kk + 1):nxxl) {
      mu[j] <- sum(para.ar * xxl[(j - 1):(j - k.ar)]) +
        sum(para.ma * epsl[(j - 1):(j - k.ma)])
    }

    # Remove added initial values
    xxl <- xxl[(kk + 1):nxxl]
    mu <- mu[(kk + 1):nxxl]

    # Calculate conditional log-likelihood
    xx2 <- (xxl - mu)^2
    l <- -0.5 * (nxxl - kk) * log(sigma^2) - 0.5 * sum(xx2 / (sigma^2))
    return(l)
  }

  # Set initial conditions for optimization
  phi.ini <- if (p == 0) NULL else rep(0.1, p)
  phiS1.ini <- if (P1 == 0) NULL else rep(0.1, P1)
  phiS2.ini <- if (P2 == 0) NULL else rep(0.1, P2)
  theta.ini <- if (q == 0) NULL else rep(-0.1, q)
  thetaS1.ini <- if (Q1 == 0) NULL else rep(-0.1, Q1)
  thetaS2.ini <- if (Q2 == 0) NULL else rep(0.1, Q2)
  sigma.ini <- sd(xx)

  beta.ini <- c(
    phi.ini,
    phiS1.ini,
    phiS2.ini,
    theta.ini,
    thetaS1.ini,
    thetaS2.ini,
    sigma.ini
  )

  # Optimize log-likelihood
  fit <- optim(
    par = beta.ini,
    fn = loglik,
    gr = NULL,
    method = "BFGS",
    control = c(list(fnscale = -1), optim_control),
    hessian = TRUE
  )

  # Check convergence
  if (fit$convergence != 0) {
    warning(paste("Optimization did not converge. Code:", fit$convergence))
  }

  # Extract estimated parameters
  npar <- length(fit$par)
  beta.hat <- c(fit$par[1:(npar - 1)], fit$par[npar]^2)
  sigma2.hat <- beta.hat[npar]

  phi.hat <- if (p != 0) beta.hat[1:p] else NULL
  phiS1.hat <- if (P1 != 0) beta.hat[(p + 1):(p + P1)] else NULL
  phiS2.hat <- if (P2 != 0) beta.hat[(p + P1 + 1):(p + P1 + P2)] else NULL
  theta.hat <- if (q != 0) {
    beta.hat[(p + P1 + P2 + 1):(p + P1 + P2 + q)]
  } else {
    NULL
  }
  thetaS1.hat <- if (Q1 != 0) {
    beta.hat[(p + P1 + P2 + q + 1):(p + P1 + P2 + q + Q1)]
  } else {
    NULL
  }
  thetaS2.hat <- if (Q2 != 0) {
    beta.hat[(p + P1 + P2 + q + Q1 + 1):(p + P1 + P2 + q + Q1 + Q2)]
  } else {
    NULL
  }

  # Calculate standard errors and p-values
  Coef <- fit$par[1:(npar - 1)]
  fit$hessian <- fit$hessian[1:(npar - 1), 1:(npar - 1)]

  # Variance-covariance matrix
  if (any(is.na(fit$hessian)) || det(-fit$hessian) <= 0) {
    warning("Hessian is singular. Standard errors may be unreliable.")
    var.coef <- matrix(NA, nrow = npar - 1, ncol = npar - 1)
    St.err <- rep(NA, npar - 1)
  } else {
    j0 <- (-fit$hessian) / nstar
    var.coef <- solve(j0) / nstar
    St.err <- sqrt(diag(var.coef))
  }

  tstat <- Coef / St.err
  pval <- 2 * pnorm(abs(tstat), lower.tail = FALSE)

  # Calculate information criteria
  npar_total <- length(Coef) + 1 # +1 for sigma^2
  loglik <- fit$value
  aic <- -2 * loglik + 2 * npar_total
  bic <- -2 * loglik + log(nstar) * npar_total
  aicc <- aic + (2 * npar_total * (npar_total + 1)) / (nstar - npar_total - 1)

  # Create coefficient names
  coef_names <- c(
    if (!is.null(phi.hat)) paste("AR", 1:length(phi.hat), sep = "."),
    if (!is.null(phiS1.hat)) paste("SAR1", 1:length(phiS1.hat), sep = "."),
    if (!is.null(phiS2.hat)) paste("SAR2", 1:length(phiS2.hat), sep = "."),
    if (!is.null(theta.hat)) paste("MA", 1:length(theta.hat), sep = "."),
    if (!is.null(thetaS1.hat)) paste("SMA1", 1:length(thetaS1.hat), sep = "."),
    if (!is.null(thetaS2.hat)) paste("SMA2", 1:length(thetaS2.hat), sep = ".")
  )

  # Create results table
  coef_table <- data.frame(
    Estimate = Coef,
    `Std. Error` = St.err,
    `z value` = tstat,
    `Pr(>|z|)` = pval,
    check.names = FALSE
  )
  rownames(coef_table) <- coef_names

  # Calculate fitted values and residuals
  # Build polynomials with estimated parameters
  poly.ar <- if (p == 0) polynomial(1) else polynomial(c(1, -phi.hat))

  poly.ar.s1 <- polynomial(1)
  if (P1 > 0) {
    for (j in 1:P1) {
      poly.ar.s1 <- polynomial(c(poly.ar.s1, rep(0, S1 - 1), -phiS1.hat[j]))
    }
  }

  poly.ar.s2 <- polynomial(1)
  if (P2 > 0) {
    for (j in 1:P2) {
      poly.ar.s2 <- polynomial(c(poly.ar.s2, rep(0, S2 - 1), -phiS2.hat[j]))
    }
  }

  poly.ma <- if (q == 0) polynomial(1) else polynomial(c(1, theta.hat))

  poly.ma.s1 <- polynomial(1)
  if (Q1 > 0) {
    for (j in 1:Q1) {
      poly.ma.s1 <- polynomial(c(poly.ma.s1, rep(0, S1 - 1), thetaS1.hat[j]))
    }
  }

  poly.ma.s2 <- polynomial(1)
  if (Q2 > 0) {
    for (j in 1:Q2) {
      poly.ma.s2 <- polynomial(c(poly.ma.s2, rep(0, S2 - 1), thetaS2.hat[j]))
    }
  }

  # Integration polynomials
  poly.i <- polynomial(c(1, -di))
  if (!is.na(S1)) {
    poly.s1.i <- polynomial(c(1, rep(0, S1 - 1), -Di1))
  } else {
    poly.s1.i <- polynomial(1)
  }

  if (!is.na(S2)) {
    poly.s2.i <- polynomial(c(1, rep(0, S2 - 1), -Di2))
  } else {
    poly.s2.i <- polynomial(1)
  }

  # Total polynomials
  poly.ar.tot <- poly.ar * poly.ar.s1 * poly.ar.s2
  poly.ma.tot <- poly.ma * poly.ma.s1 * poly.ma.s2
  poly <- poly.ar * poly.ar.s1 * poly.ar.s2 * poly.i * poly.s1.i * poly.s2.i

  # Extract parameters
  para.ar <- if ((p + P1 + P2) == 0) 0 else -poly.ar.tot[2:length(poly.ar.tot)]
  para.ma <- if ((q + Q1 + Q2) == 0) 0 else poly.ma.tot[2:length(poly.ma.tot)]
  para <- if ((p + P1 + P2 + di + Di1 + Di2) == 0) 0 else -poly[2:length(poly)]

  k.ar <- length(para.ar)
  k.ma <- length(para.ma)
  k <- length(para)

  # Prepare extended series
  nxx <- length(xx)
  kk <- max(k.ar, k.ma)
  epsl <- rep(0, nxx + kk)
  xxl <- c(rep(mean(xx[1:10]), kk), xx)
  nxxl <- length(xxl)

  # Calculate residuals on differenced series
  for (j in (kk + 1):nxxl) {
    epsl[j] <- xxl[j] -
      (sum(para.ar * xxl[(j - 1):(j - k.ar)]) +
        sum(para.ma * epsl[(j - 1):(j - k.ma)]))
  }

  # Calculate conditional means
  mu <- rep(0, nxxl)
  for (j in (kk + 1):nxxl) {
    mu[j] <- sum(para.ar * xxl[(j - 1):(j - k.ar)]) +
      sum(para.ma * epsl[(j - 1):(j - k.ma)])
  }

  # Remove added initial values
  xxl <- xxl[(kk + 1):nxxl]
  mu <- mu[(kk + 1):nxxl]

  # Calculate residuals as difference between differenced series and conditional mean
  resid <- xxl - mu
  numd <- sum(c(di, Di1 * S1, Di2 * S2), na.rm = TRUE)
  resid <- c(rep(0, numd), resid)

  fitted <- origdata - resid

  # Create output object with class
  out <- structure(
    list(
      # Model specification
      order = c(p, di, q),
      seasonal = list(order = c(P1, Di1, Q1), period = S1),
      seasonal2 = list(order = c(P2, Di2, Q2), period = S2),

      # Coefficients
      coef = setNames(Coef, coef_names),
      var.coef = var.coef,

      # Model fit statistics
      loglik = loglik,
      aic = aic,
      bic = bic,
      aicc = aicc,
      sigma2 = sigma2.hat,

      # Series and residuals
      x = origdata,
      fitted = fitted,
      residuals = resid,

      # Additional info
      series = deparse(substitute(y)),
      nobs = n,
      n.used = nstar,
      call = match.call(),
      convergence = fit$convergence
    ),
    class = "msarima"
  )

  return(out)
}

#' @rdname msarima
#' @method print msarima
#' @param x An object of class \code{"msarima"}.
#' @param ... Additional arguments (currently not used).
#' @export
print.msarima <- function(x, ...) {
  cat("Series:", x$series, "\n")
  cat(sprintf("MSARIMA(%d,%d,%d)", x$order[1], x$order[2], x$order[3]))

  if (!is.na(x$seasonal$period)) {
    cat(sprintf(
      "(%d,%d,%d)[%d]",
      x$seasonal$order[1],
      x$seasonal$order[2],
      x$seasonal$order[3],
      x$seasonal$period
    ))
  }

  if (!is.na(x$seasonal2$period)) {
    cat(sprintf(
      "(%d,%d,%d)[%d]",
      x$seasonal2$order[1],
      x$seasonal2$order[2],
      x$seasonal2$order[3],
      x$seasonal2$period
    ))
  }

  cat("\n\nCoefficients:\n")
  print(round(x$coef, 4))

  cat(sprintf("\nsigma^2 = %.4f", x$sigma2))
  cat(sprintf("  log likelihood = %.2f", x$loglik))
  cat(sprintf(
    "\nAIC = %.2f   AICc = %.2f   BIC = %.2f\n",
    x$aic,
    x$aicc,
    x$bic
  ))

  if (x$convergence != 0) {
    cat(sprintf(
      "\nWarning: Optimization convergence code = %d\n",
      x$convergence
    ))
  }

  invisible(x)
}


#' @rdname msarima
#' @method summary msarima
#' @param object An object of class \code{"msarima"}.
#' @export
summary.msarima <- function(object, ...) {
  # Calculate coefficient table
  se <- sqrt(diag(object$var.coef))
  tstat <- object$coef / se
  pval <- 2 * pnorm(abs(tstat), lower.tail = FALSE)

  coef_matrix <- cbind(
    Estimate = object$coef,
    `Std. Error` = se,
    `z value` = tstat,
    `Pr(>|z|)` = pval
  )

  out <- structure(
    list(
      call = object$call,
      series = object$series,
      order = object$order,
      seasonal = object$seasonal,
      seasonal2 = object$seasonal2,
      coefficients = coef_matrix,
      sigma2 = object$sigma2,
      loglik = object$loglik,
      aic = object$aic,
      bic = object$bic,
      aicc = object$aicc,
      nobs = object$nobs,
      n.used = object$n.used,
      convergence = object$convergence
    ),
    class = "summary.msarima"
  )

  return(out)
}

#' @rdname msarima
#' @method print summary.msarima
#' @param digits Number of digits to display. Default is 4.
#' @param signif.stars Logical; if \code{TRUE}, significance stars are printed. Default is \code{TRUE}.
#' @export
print.summary.msarima <- function(x, digits = 4, signif.stars = TRUE, ...) {
  cat("Series:", x$series, "\n")
  cat(sprintf("MSARIMA(%d,%d,%d)", x$order[1], x$order[2], x$order[3]))

  if (!is.na(x$seasonal$period)) {
    cat(sprintf(
      "(%d,%d,%d)[%d]",
      x$seasonal$order[1],
      x$seasonal$order[2],
      x$seasonal$order[3],
      x$seasonal$period
    ))
  }

  if (!is.na(x$seasonal2$period)) {
    cat(sprintf(
      "(%d,%d,%d)[%d]",
      x$seasonal2$order[1],
      x$seasonal2$order[2],
      x$seasonal2$order[3],
      x$seasonal2$period
    ))
  }

  cat("\n\nCoefficients:\n")
  printCoefmat(
    x$coefficients,
    digits = digits,
    signif.stars = signif.stars,
    P.values = TRUE,
    has.Pvalue = TRUE
  )

  cat(sprintf("\nsigma^2 = %.*f", digits, x$sigma2))
  cat(sprintf("  log likelihood = %.*f", digits, x$loglik))
  cat(sprintf(
    "\nAIC = %.*f   AICc = %.*f   BIC = %.*f",
    digits,
    x$aic,
    digits,
    x$aicc,
    digits,
    x$bic
  ))
  cat(sprintf("\nn.obs = %d   n.used = %d\n", x$nobs, x$n.used))

  if (x$convergence != 0) {
    cat(sprintf(
      "\nWarning: Optimization convergence code = %d\n",
      x$convergence
    ))
  }

  invisible(x)
}

#' @rdname msarima
#' @method coef msarima
#' @export
coef.msarima <- function(object, ...) {
  return(object$coef)
}

#' @rdname msarima
#' @method fitted msarima
#' @export
fitted.msarima <- function(object, ...) {
  return(object$fitted)
}

#' @rdname msarima
#' @method residuals msarima
#' @export
residuals.msarima <- function(object, ...) {
  return(object$residuals)
}

#' @rdname msarima
#' @method logLik msarima
#' @export
logLik.msarima <- function(object, ...) {
  structure(
    object$loglik,
    df = length(object$coef) + 1,
    nobs = object$n.used,
    class = "logLik"
  )
}

#' @rdname msarima
#' @method AIC msarima
#' @export
AIC.msarima <- function(object, ...) {
  return(object$aic)
}

#' @rdname msarima
#' @method BIC msarima
#' @export
BIC.msarima <- function(object, ...) {
  return(object$bic)
}


check_seasonal <- function(x, name = deparse(substitute(x))) {
  # Object must be a list
  if (!is.list(x)) {
    stop(sprintf("'%s' must be a list", name))
  }

  # List must contain exactly two elements: order and period
  required_names <- c("order", "period")
  if (!all(required_names %in% names(x)) || length(x) != 2) {
    stop(sprintf(
      "'%s' must be a list with exactly two elements: list(order = c(0, 0, 0), period = NA)",
      name
    ))
  }

  # Check 'order'
  if (!is.numeric(x$order) || length(x$order) != 3) {
    stop(sprintf("'%s$order' must be a numeric vector of length 3", name))
  }
  if (any(is.na(x$order))) {
    stop(sprintf("'%s$order' must not contain NA values", name))
  }
  if (any(x$order < 0) || any(x$order != floor(x$order))) {
    stop(sprintf("'%s$order' must contain non-negative integers", name))
  }

  # Check 'period'
  if (is.na(x$period)) {
    # Constraint: if period is NA, then order must be c(0, 0, 0)
    if (any(x$order != 0)) {
      stop(sprintf(
        "If '%s$period' is NA then '%s$order' must be c(0, 0, 0)",
        name,
        name
      ))
    }
  } else {
    # Period must be a single positive integer
    if (!is.numeric(x$period) || length(x$period) != 1 || is.na(x$period)) {
      stop(sprintf("'%s$period' must be NA or a single numeric value", name))
    }
    if (x$period <= 0 || x$period != floor(x$period)) {
      stop(sprintf("'%s$period' must be a positive integer", name))
    }
  }

  invisible(TRUE)
}
