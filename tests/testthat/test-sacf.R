test_that("sacf picks up the seasonal periods and highlights the right lags", {
  y <- msts(rnorm(500), seasonal.periods = c(7, 30))

  cc <- sacf(y, plot = FALSE)
  expect_equal(attr(cc, "season"), c(7, 30))
  expect_equal(max(cc$lag), 90) # 3 * max(season)
  expect_equal(drop(cc$acf), drop(acf(as.numeric(y), 90, plot = FALSE)$acf))

  grp <- attr(cc, "season.group")
  expect_equal(grp[cc$lag %in% c(0, 7, 30, 60)], c(NA, 1L, 2L, 2L))

  # a lag that is a multiple of both periods goes to the longest one
  cc2 <- sacf(y, season = c(7, 14), plot = FALSE)
  grp2 <- attr(cc2, "season.group")
  expect_equal(grp2[cc2$lag %in% c(7, 14, 21, 28)], c(1L, 2L, 1L, 2L))
})

test_that("sacf/spacf on a fitted model drop the leading NA residuals", {
  y <- msarima.sim(nsim = 300, order = c(1, 0, 0), seasonal = list(c(0, 1, 1)),
                   seasonal.period = 7, ar = 0.5, ma = 0.4, sigma2 = 1)
  fit <- msarima(y, order = c(1, 0, 0), seasonal = list(c(0, 1, 1)),
                 seasonal.period = 7)

  cc <- spacf(fit, plot = FALSE)
  expect_equal(attr(cc, "season"), 7)
  expect_equal(cc$n.used, fit$nobs - 7)
  expect_false(anyNA(cc$acf))
})

test_that("no seasonal period is not an error", {
  cc <- sacf(rnorm(200), season = numeric(0), plot = FALSE)
  expect_length(attr(cc, "season"), 0)
})

test_that("ggsacf draws the non-seasonal lags too", {
  skip_if_not_installed("ggplot2")
  y <- msts(rnorm(500), seasonal.periods = c(7, 30))

  b <- ggplot2::ggplot_build(ggsacf(y, lag.max = 60))
  n <- sapply(b$data, nrow)
  expect_true(any(n == 61)) # all the lags, 0 to 60
  expect_true(any(n == 10)) # the seasonal ones: 7, 14, 21, 28, 30, 35, ...
})
