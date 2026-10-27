#' Multi-Seasonal Portmanteau Tests
#'
#' @description
#' Multi-seasonal generalisation of the Ljung-Box and Monti portmanteau tests:
#' residual autocorrelation is checked at low-order and seasonal lags at once,
#' for an arbitrary number of seasonal periods.
#'
#' @param x A numeric vector, a `ts`/`msts` object, or a fitted [msarima] model
#'   (in which case the model residuals are used).
#' @param type Either `"Ljung-Box"` (default), which builds the statistic from
#'   the autocorrelations, or `"Monti"`, which uses the partial
#'   autocorrelations. Partial matching is allowed.
#' @param season Numeric vector of seasonal periods to test, of any length. If
#'   `NULL` (default) the periods are taken from the `msts`/`ts` attributes of
#'   `x`, or from the fitted model.
#' @param nlag Numeric vector giving the number of seasonal multiples of each
#'   period to test, i.e. lags `season * (1:nlag)`. A single value is recycled across
#'   `season`; otherwise give one value per period, e.g. `c(5, 3)`. Default `5`.
#' @param lags_short Number of low-order, non-seasonal lags `1:lags_short` added
#'   to the seasonal ones. Default `0` (none).
#' @param npar Number of estimated model parameters. Default `0` for a raw
#'   series; for a fitted [msarima] model it defaults to the number of estimated
#'   (non-fixed) AR, seasonal AR, MA and seasonal MA coefficients.
#'
#' @details
#' The lag set tested is `sort(unique(c(1:lags_short, outer(season, 1:nlag))))`,
#' the same construction as in the original two-season versions of these
#' tests, generalised to any number of periods.
#'
#' As in [sacf], leading and trailing `NA`s (such as those produced by
#' differencing in [msarima]) are dropped ([stats::na.contiguous]) before
#' computing the statistic, and the number of observations used is the length
#' of the series that remains.
#'
#' @return An object of class `"htest"`, with:
#' \item{statistic}{the test statistic ("X-squared").}
#' \item{parameter}{the degrees of freedom ("df"), the number of lags tested
#'   minus `npar`.}
#' \item{p.value}{the p-value, from the upper tail of a chi-squared
#'   distribution with `parameter` degrees of freedom.}
#' \item{method}{a character string describing the test.}
#' \item{data.name}{a character string naming the data.}
#' \item{lag}{the lags actually tested.}
#'
#' @examples
#' fit <- msarima(MAD_m40, order = c(1, 0, 2),
#'                seasonal = list(c(0, 1, 1), c(0, 1, 1)))
#'
#' # periods 24 and 168 taken from the msts object
#' msportmanteau(MAD_m40, nlag = 3)
#' # residual test, periods and npar from the model
#' msportmanteau(fit, type = "Monti", nlag = 3)
#'
#' @references
#' Grigoletto, M. and Lisi, F. (2026). Portmanteau goodness-of-fit tests in the
#' presence of simple, mixed, and multiple seasonal autocorrelation. *arXiv*.
#' \doi{10.48550/arXiv.2609.32381}
#'
#' Ljung, G. M. and Box, G. E. P. (1978). On a measure of lack of fit in time
#' series models. *Biometrika*, 65(2), 297-303.
#'
#' Monti, A. C. (1994). A proposal for a residual autocorrelation test in
#' linear models. *Biometrika*, 81(4), 776-780.
#'
#' @seealso [sacf], [msarima], [stats::Box.test]
#'
#' @importFrom stats acf pchisq na.contiguous
#'
#' @export
msportmanteau <- function(
  x,
  type = c("Ljung-Box", "Monti"),
  season = NULL,
  nlag = 5,
  lags_short = 0,
  npar = NULL
) {
  data.name <- deparse1(substitute(x))
  type <- match.arg(type)

  if (inherits(x, "msarima")) {
    data.name <- x$series
    if (is.null(season)) {
      season <- x$seasonal.period
    }
    if (is.null(npar)) {
      npar <- sum(!x$mask[names(x$coef) != "intercept"])
    }
    x <- x$residuals
  } else if (is.null(season)) {
    season <- if (inherits(x, "msts")) {
      attr(x, "msts")
    } else if (is.ts(x)) {
      frequency(x)
    }
  }
  if (is.null(npar)) {
    npar <- 0
  }

  season <- sort(unique(season[!is.na(season) & season > 1]))
  if (length(season) == 0) {
    stop("'season' must contain at least one seasonal period greater than 1")
  }

  nlag <- rep_len(nlag, length(season))
  lags_seasonal <- unlist(mapply(
    function(s, m) s * seq_len(m),
    season,
    nlag,
    SIMPLIFY = FALSE
  ))
  lags_short <- if (lags_short > 0) seq_len(lags_short) else numeric(0)
  L <- sort(unique(c(lags_short, lags_seasonal)))

  xc <- stats::na.contiguous(as.numeric(x))
  n <- length(xc)

  if (max(L) > n - 1) {
    stop(sprintf(
      "The requested lags go up to %d, but only %d observations are available after removing leading/trailing NAs: the maximum usable lag is %d",
      max(L),
      n,
      n - 1
    ))
  }

  PARAMETER <- length(L) - npar
  if (PARAMETER <= 0) {
    stop(sprintf(
      "'npar' (%d) leaves no degrees of freedom with %d lags tested",
      npar,
      length(L)
    ))
  }

  partial <- type == "Monti"
  cc <- stats::acf(
    xc,
    lag.max = max(L),
    type = if (partial) "partial" else "correlation",
    plot = FALSE
  )
  # the acf starts at lag 0, the pacf at lag 1
  vals <- drop(cc$acf)^2
  vals <- vals[if (partial) L else L + 1]

  STATISTIC <- n * (n + 2) * sum(vals / (n - L))
  PVALUE <- stats::pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)

  structure(
    list(
      statistic = c("X-squared" = STATISTIC),
      parameter = c(df = PARAMETER),
      p.value = PVALUE,
      method = paste("Multi-seasonal", type, "test"),
      data.name = data.name,
      lag = L
    ),
    class = "htest"
  )
}
