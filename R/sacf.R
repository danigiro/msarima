#' Autocorrelation Plots for Multiple Seasonal Cycles
#'
#' @description
#' Sample autocorrelation (`sacf`) and partial autocorrelation (`spacf`) functions
#' that highlight the lags belonging to one or more seasonal cycles, so that the
#' seasonal structure of a multi-seasonal series is visible at a glance.
#' `ggsacf()` and `ggspacf()` are the \pkg{ggplot2} equivalents.
#'
#' @param x A numeric vector, a `ts`/`msts` object, or a fitted [msarima] model
#'   (in which case the model residuals are used).
#' @param season Numeric vector of seasonal periods to highlight, of any length.
#'   If `NULL` (default) the periods are taken from the `msts`/`ts` attributes of
#'   `x`, or from the fitted model.
#' @param lag.max Maximum lag. Defaults to `3 * max(season)`, or to the
#'   [stats::acf] default when no seasonal period is available.
#' @param ci Coverage of the confidence bands, as a probability. Default `0.95`.
#' @param col Colours used for the seasonal lags, one per element of `season`
#'   (recycled). Defaults to a colour-blind safe palette (Okabe and Ito, 2008),
#'   the same one used by `ggsacf()` and `ggspacf()`.
#' @param main,ylab Plot title and y-axis label.
#' @param plot If `TRUE` (default) the correlogram is drawn.
#' @param ... Further arguments passed to [stats::acf], and from `spacf()` and
#'   `ggspacf()` to `sacf()` and `ggsacf()`.
#'
#' @details
#' Lags that are multiples of more than one period (lag 168 with
#' `season = c(24, 168)`, say) are attributed to the longest matching period.
#'
#' The series is treated as a plain numeric vector, so lags are always reported in
#' number of observations and not in units of time. Leading and trailing `NA`s,
#' such as those produced by differencing in [msarima], are dropped
#' ([stats::na.contiguous]).
#'
#' @return `sacf()` and `spacf()` invisibly return the [stats::acf] object, with
#'   the seasonal periods in the `"season"` attribute and, for each lag, the index
#'   of the period it belongs to in `"season.group"`. `ggsacf()` and `ggspacf()`
#'   return a `ggplot` object.
#'
#' @examples
#' fit <- msarima(MAD_m40, order = c(1, 0, 2),
#'                seasonal = list(c(0, 1, 1), c(0, 1, 1)))
#'
#' sacf(MAD_m40, lag.max = 200)   # periods 24 and 168 taken from the msts object
#' spacf(fit)                     # residual pacf, periods taken from the model
#'
#' @seealso [stats::acf], [msarima]
#'
#' @name sacf
#' @export
sacf <- function(x, season = NULL, lag.max = NULL, ci = 0.95, col = NULL,
                 main = NULL, ylab = "ACF", plot = TRUE, ...) {
  if (inherits(x, "msarima")) {
    if (is.null(season)) season <- x$seasonal.period
    x <- x$residuals
  } else if (is.null(season)) {
    season <- if (inherits(x, "msts")) attr(x, "msts") else if (is.ts(x)) frequency(x)
  }
  season <- sort(unique(season[!is.na(season) & season > 1]))

  x <- as.numeric(x)
  if (is.null(lag.max)) {
    lag.max <- if (length(season) > 0) 3 * max(season) else 10 * log10(length(x))
  }
  cc <- stats::acf(x, lag.max = lag.max, plot = FALSE,
                   na.action = stats::na.contiguous, ...)

  # For each lag, the period it belongs to: the longest one, since `season` is sorted
  lag <- drop(cc$lag)
  group <- rep(NA_integer_, length(lag))
  for (i in seq_along(season)) group[lag > 0 & lag %% season[i] == 0] <- i
  attr(cc, "season") <- season
  attr(cc, "season.group") <- group

  if (plot) {
    val <- drop(cc$acf)
    clim <- qnorm((1 + ci) / 2) / sqrt(cc$n.used)
    plot(lag, val, type = "h", lwd = 1.2, xlab = "Lag", ylab = ylab, main = main,
         ylim = range(val, clim, -clim))
    abline(h = 0)
    abline(h = c(-clim, clim), col = "grey40", lty = 2)
    if (is.null(col)) col <- c("#E69F00", "#0072B2", "#009E73", "#D55E00", "#CC79A7")
    col <- rep_len(col, max(length(season), 1L))
    for (i in seq_along(season)) {
      k <- which(group == i)
      points(lag[k], val[k], type = "h", lwd = 1.2, col = col[i])
      points(lag[k], val[k], pch = 20, cex = 1.1, col = col[i])
    }
  }
  invisible(cc)
}

#' @rdname sacf
#' @export
spacf <- function(x, ..., ylab = "PACF") {
  sacf(x, type = "partial", ylab = ylab, ...)
}

#' @rdname sacf
#' @export
ggsacf <- function(x, season = NULL, lag.max = NULL, ci = 0.95,
                   main = NULL, ylab = "ACF", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function")
  }
  cc <- sacf(x, season = season, lag.max = lag.max, plot = FALSE, ...)
  season <- attr(cc, "season")
  clim <- qnorm((1 + ci) / 2) / sqrt(cc$n.used)
  df <- data.frame(
    lag = drop(cc$lag),
    acf = drop(cc$acf),
    season = factor(season[attr(cc, "season.group")], levels = season)
  )
  # All the lags first, then the seasonal ones on top, as in sacf()
  dfs <- df[!is.na(df$season), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(yintercept = c(-clim, clim), colour = "grey40", linetype = 2) +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::geom_segment(data = dfs, ggplot2::aes(xend = lag, yend = 0, colour = season)) +
    ggplot2::geom_point(data = dfs, ggplot2::aes(colour = season), size = 1.5) +
    ggplot2::scale_colour_manual(
      values = c("#E69F00", "#0072B2", "#009E73", "#D55E00", "#CC79A7")) +
    ggplot2::labs(x = "Lag", y = ylab, title = main, colour = "Period")
}

#' @rdname sacf
#' @export
ggspacf <- function(x, ..., ylab = "PACF") {
  ggsacf(x, type = "partial", ylab = ylab, ...)
}
