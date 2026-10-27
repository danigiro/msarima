#' Hourly Traffic Volume on the M-40 Ring Road, Madrid
#'
#' @description
#' Number of vehicles per hour recorded by a single traffic sensor on the M-40, the
#' outer ring road of Madrid, over the eight weeks from 3 February to 29 March 2024.
#'
#' The series shows the two seasonal patterns that \code{\link{msarima}} is meant for:
#' a daily cycle of period 24, with the morning and evening commuting peaks, and a
#' weekly cycle of period 168, with markedly lower volumes at the weekend.
#'
#' @format
#' An \code{msts} object of length 1344 (56 days, or eight complete weeks, of hourly
#' data) with seasonal periods 24 and 168. The first observation is midnight of
#' Saturday 3 February 2024 and the last is 11 p.m. of Friday 29 March 2024.
#'
#' @details
#' The raw measurements are half-hourly counts from sensor 3753, recorded between
#' 1 February and 31 May 2024 and summed in pairs to obtain hourly volumes; the eight
#' weeks kept here are the window studied by Grigoletto and Lisi (2026), who chose this
#' sensor for having the fewest missing observations in the network. A few of the
#' underlying half-hourly counts had been interpolated by the provider, so 3 of the
#' hourly values are not whole numbers. There are no missing values and no gaps in the
#' hourly grid.
#'
#' The script that builds the object from the raw file is available at
#' \url{https://github.com/danigiro/msarima/blob/main/data-raw/MAD_m40.R}.
#'
#' @references
#' Grigoletto, M. and Lisi, F. (2026). Portmanteau goodness-of-fit tests in the
#' presence of simple, mixed, and multiple seasonal autocorrelation. *arXiv*.
#' \doi{10.48550/arXiv.2609.32381}
#'
#' @source
#' Madrid City Council open data portal, \url{https://datos.madrid.es}.
#'
#' @examples
#' # The daily and weekly cycles
#' plot(window(MAD_m40, end = 4), ylab = "vehicles per hour")
#'
#' # A multi-seasonal ARIMA with one daily and one weekly difference
#' \donttest{
#' fit <- msarima(MAD_m40,
#'                order = c(1, 0, 2),
#'                seasonal = list(c(0, 1, 1), c(0, 1, 1)))
#' summary(fit)
#' plot(forecast(fit, h = 48))
#' }
"MAD_m40"
