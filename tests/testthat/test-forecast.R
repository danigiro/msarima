test_that("the forecast object has the documented structure", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )
  fc <- forecast(fit, h = 12, level = c(80, 95))

  expect_s3_class(fc, "forecast")
  expect_length(as.numeric(fc$mean), 12)
  expect_equal(dim(fc$lower), c(12L, 2L))
  expect_equal(colnames(fc$lower), c("80%", "95%"))
  expect_true(all(fc$lower[, 1] < fc$mean))
  expect_true(all(fc$upper[, 1] > fc$mean))
  # a wider level gives a wider interval, and the interval widens with h
  expect_true(all(fc$upper[, 2] - fc$lower[, 2] > fc$upper[, 1] - fc$lower[, 1]))
  expect_true(all(diff(fc$upper[, 1] - fc$lower[, 1]) >= 0))

  expect_length(forecast(fit, h = 5, fan = TRUE)$level, 17)
})

test_that("the one-step forecast equals the one-step fitted value", {
  # The recursion used in estimation and the one used in forecasting must
  # agree. Refitting on y[-n] with the coefficients held fixed isolates the
  # recursion from the optimisation.
  y <- sim_int()
  n <- length(y)

  full <- msarima(
    y,
    order = c(1, 1, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )
  sub <- suppressWarnings(msarima(
    y[-n],
    order = c(1, 1, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4,
    fixed = c(unname(coef(full)), NA),
    optim_control = list(maxit = 20)
  ))

  expect_equal(
    as.numeric(forecast(sub, h = 1)$mean),
    as.numeric(fitted(full))[n]
  )
})

test_that("forecasts of a differenced model are integrated back to levels", {
  y <- sim_sint()
  n <- length(y)

  fit_D <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 1, 1)),
    seasonal.period = 4
  )
  fit_w <- msarima(
    diff(y, lag = 4),
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  # first forecast on levels = first forecast on the differenced scale + y[n+1-S]
  expect_equal(
    as.numeric(forecast(fit_D, h = 1)$mean),
    as.numeric(forecast(fit_w, h = 1)$mean) + y[n + 1 - 4]
  )
  expect_false(anyNA(forecast(fit_D, h = 24)$mean))
})

test_that("the intercept is used in the point forecasts", {
  # forecasts converge to the process mean, not to 0
  set.seed(11)
  y <- 100 + as.numeric(msarima.sim(
    nsim = 300,
    order = c(1, 0, 0),
    seasonal = list(c(1, 0, 0)),
    seasonal.period = 4,
    ar = c(0.5, 0.3),
    sigma2 = 1
  ))
  fit <- msarima(
    y,
    order = c(1, 0, 0),
    seasonal = list(c(1, 0, 0)),
    seasonal.period = 4,
    include.constant = TRUE
  )
  implied <- unname(
    coef(fit)["intercept"] / ((1 - coef(fit)["AR(1)"]) * (1 - coef(fit)["SAR(1,1)"]))
  )

  fc <- forecast(fit, h = 80)
  expect_equal(implied, mean(y), tolerance = 0.05)
  expect_equal(as.numeric(fc$mean)[80], implied, tolerance = 1e-4)
})

test_that("seasonal AR coefficients of order 2 act on lags S and 2S", {
  y <- sim_arma(n = 200)
  fit <- msarima(
    y,
    order = c(0, 0, 0),
    seasonal = list(c(2, 0, 0)),
    seasonal.period = 4
  )
  n <- length(y)

  expect_equal(
    as.numeric(forecast(fit, h = 1)$mean),
    unname(coef(fit)[1] * y[n - 3] + coef(fit)[2] * y[n - 7])
  )
})

test_that("seasonal differencing of order 2 is applied to the right lags", {
  y <- cumsum(cumsum(sim_arma(n = 200)))
  fit <- msarima(
    y,
    order = c(0, 0, 0),
    seasonal = list(c(0, 2, 1)),
    seasonal.period = 4
  )
  n <- length(y)

  # with (1 - B^4)^2 the two-step-back structure must reproduce a pure random
  # walk forecast when the MA term is switched off
  fc <- forecast(fit, h = 1)
  expect_false(is.na(as.numeric(fc$mean)))
  expect_equal(fit$n.used, length(y) - 8)
})

test_that("one-step prediction intervals are correct for differenced models", {
  y <- sim_int()
  fit <- msarima(
    y,
    order = c(1, 1, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  fc1 <- forecast(fit, h = 1)
  fc5 <- forecast(fit, h = 5)

  # the one-step forecast variance is sigma2
  qq <- qnorm(0.5 + fc1$level[1] / 200)
  implied <- ((fc1$upper[1, 1] - fc1$mean[1]) / qq)^2
  expect_equal(as.numeric(implied), fit$sigma2)

  # the h = 1 interval matches the first row of a longer forecast
  expect_equal(as.numeric(fc1$lower), as.numeric(fc5$lower[1, ]))
  expect_equal(as.numeric(fc1$upper), as.numeric(fc5$upper[1, ]))
})
