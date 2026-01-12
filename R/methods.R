#' Methods for msarima Objects
#'
#' @description
#' Methods for extracting information from and displaying fitted multi-seasonal ARIMA models.
#'
#' @param x,object An object of class \code{"msarima"}.
#' @param digits Number of digits to display. Default is 4.
#' @param signif.stars Logical; if \code{TRUE}, significance stars are printed. Default is \code{TRUE}.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' \describe{
#'   \item{\code{print}}{Prints a concise summary of the fitted model including coefficients
#'     and information criteria.}
#'   \item{\code{summary}}{Produces a detailed summary including coefficient estimates,
#'     standard errors, z-values, and p-values.}
#'   \item{\code{coef}}{Extracts the estimated model coefficients.}
#'   \item{\code{fitted}}{Extracts the fitted values.}
#'   \item{\code{residuals}}{Extracts the model residuals.}
#'   \item{\code{logLik}}{Extracts the log-likelihood value.}
#'   \item{\code{AIC}}{Calculates the Akaike Information Criterion.}
#'   \item{\code{BIC}}{Calculates the Bayesian Information Criterion.}
#' }
#'
#' @return
#' \describe{
#'   \item{\code{print}}{Returns the input object invisibly.}
#'   \item{\code{summary}}{Returns an object of class \code{"summary.msarima"} containing:
#'     \itemize{
#'       \item \code{call}: The matched call.
#'       \item \code{series}: Name of the time series.
#'       \item \code{order}: Non-seasonal ARIMA order.
#'       \item \code{seasonal}: List of seasonal component specifications.
#'       \item \code{seasonal.period}: Vector of seasonal periods.
#'       \item \code{coefficients}: Matrix with estimates, standard errors, z-values, and p-values.
#'       \item \code{sigma2}: Estimated variance of the residuals.
#'       \item \code{loglik}: Log-likelihood of the fitted model.
#'       \item \code{aic}: Akaike Information Criterion.
#'       \item \code{bic}: Bayesian Information Criterion.
#'       \item \code{aicc}: Corrected Akaike Information Criterion.
#'       \item \code{nobs}: Number of observations in the original series.
#'       \item \code{n.used}: Number of observations used in estimation.
#'       \item \code{method}: Description of the fitted model.
#'       \item \code{convergence}: Convergence code from optimization.
#'     }
#'   }
#'   \item{\code{coef}}{A named numeric vector of coefficients.}
#'   \item{\code{fitted}}{A time series object of fitted values.}
#'   \item{\code{residuals}}{A time series object of residuals.}
#'   \item{\code{logLik}}{An object of class \code{"logLik"}.}
#'   \item{\code{AIC}}{The AIC value.}
#'   \item{\code{BIC}}{The BIC value.}
#' }
#'
#' @seealso \code{\link{msarima}}
#'
#' @examples
#' # Fit a model
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
#' # Use methods
#' print(fit)
#' summary(fit)
#' coef(fit)
#' head(fitted(fit))
#' head(residuals(fit))
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#'
#' @name msarima-methods
NULL

#' @rdname msarima-methods
#' @method print msarima
#' @export
print.msarima <- function(x, ...) {
  cat("Series:", x$series, "\n")
  cat(x$method)

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


#' @rdname msarima-methods
#' @method summary msarima
#' @export
summary.msarima <- function(object, ...) {
  # Calculate coefficient table
  se <- sqrt(diag(object$var_coef))
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
      seasonal.period = object$seasonal.period,
      coefficients = coef_matrix,
      sigma2 = object$sigma2,
      loglik = object$loglik,
      aic = object$aic,
      bic = object$bic,
      aicc = object$aicc,
      nobs = object$nobs,
      n.used = object$n.used,
      method = object$method,
      convergence = object$convergence
    ),
    class = "summary.msarima"
  )

  return(out)
}

#' @rdname msarima-methods
#' @method print summary.msarima
#' @export
print.summary.msarima <- function(x, digits = 4, signif.stars = TRUE, ...) {
  cat("Series:", x$series, "\n")
  cat(x$method)

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

  if (x$convergence != 0) {
    cat(sprintf(
      "\nWarning: Optimization convergence code = %d\n",
      x$convergence
    ))
  }

  invisible(x)
}

#' @rdname msarima-methods
#' @method coef msarima
#' @export
coef.msarima <- function(object, ...) {
  return(object$coef)
}

#' @rdname msarima-methods
#' @method fitted msarima
#' @export
fitted.msarima <- function(object, ...) {
  return(object$fitted)
}

#' @rdname msarima-methods
#' @method residuals msarima
#' @export
residuals.msarima <- function(object, ...) {
  return(object$residuals)
}

#' @rdname msarima-methods
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

#' @rdname msarima-methods
#' @method AIC msarima
#' @export
AIC.msarima <- function(object, ...) {
  return(object$aic)
}

#' @rdname msarima-methods
#' @method BIC msarima
#' @export
BIC.msarima <- function(object, ...) {
  return(object$bic)
}
