#' Lag polynomial 1 + sign * (b_1 L^s + ... + b_p L^{p s}), multiplied over all components
#'
#' @param coefs Coefficients in the user's order.
#' @param orders Order per component (non-seasonal first).
#' @param periods \code{c(1, seasonal.period)}.
#' @param sign -1 for AR, +1 for MA.
#' @noRd
build_lag_poly <- function(coefs, orders, periods, sign = 1) {
  if (sum(orders) == 0) return(polynomial(1))
  groups <- split(coefs, rep(seq_along(orders), orders))
  do.call(prod, Map(function(b, s) {
    p <- numeric(length(b) * s)
    p[seq_along(b) * s] <- b
    polynomial(c(1, sign * p))
  }, groups, periods[orders > 0]))
}

#' Rebuild y as a ts/msts series with the frequency and periods of x
#' @noRd
ts_like <- function(x, y, start) {
  if (inherits(x, "msts")) {
    msts(y, seasonal.periods = attr(x, "msts"), ts.frequency = frequency(x),
         start = start)
  } else {
    ts(y, start = start, frequency = frequency(x))
  }
}

#' y aligned with x, padded with leading NAs (fitted values, residuals)
#' @noRd
align_ts <- function(x, y) {
  stopifnot(length(y) <= length(x))
  ts_like(x, c(rep(NA, length(x) - length(y)), y), start(x))
}

#' y starting right after the end of x (forecasts)
#' @noRd
ahead_ts <- function(x, y) {
  ts_like(x, y, tsp(x)[2] + 1 / frequency(x))
}
