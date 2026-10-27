test_that("d = 1 is equivalent to fitting the ARMA on diff(y)", {
  y <- sim_int()

  fit_d <- msarima(
    y,
    order = c(1, 1, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )
  fit_w <- msarima(
    diff(y),
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  expect_equal(unname(coef(fit_d)), unname(coef(fit_w)))
  expect_equal(fit_d$sigma2, fit_w$sigma2)
  expect_equal(fit_d$loglik, fit_w$loglik)
  expect_equal(
    as.numeric(residuals(fit_d))[-1],
    as.numeric(residuals(fit_w))
  )
})

test_that("seasonal differencing is equivalent to fitting on diff(y, S)", {
  y <- sim_sint()

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

  expect_equal(unname(coef(fit_D)), unname(coef(fit_w)))
  expect_equal(fit_D$sigma2, fit_w$sigma2)
  expect_equal(
    as.numeric(residuals(fit_D))[-(1:4)],
    as.numeric(residuals(fit_w))
  )
})

test_that("non-seasonal and seasonal differencing combine correctly", {
  y <- cumsum(sim_sint())

  fit_dD <- msarima(
    y,
    order = c(1, 1, 1),
    seasonal = list(c(0, 1, 1)),
    seasonal.period = 4
  )
  fit_w <- msarima(
    diff(diff(y), lag = 4),
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  expect_equal(unname(coef(fit_dD)), unname(coef(fit_w)))
  expect_equal(
    as.numeric(residuals(fit_dD))[-(1:5)],
    as.numeric(residuals(fit_w))
  )
})

test_that("residuals and fitted are NA-padded and aligned with the input", {
  y <- sim_int()
  fit <- msarima(
    y,
    order = c(1, 1, 1),
    seasonal = list(c(0, 1, 1)),
    seasonal.period = 4
  )
  n_diff <- 1 + 4 # d + sum(D * S)

  expect_length(as.numeric(residuals(fit)), length(y))
  expect_length(as.numeric(fitted(fit)), length(y))
  expect_equal(sum(is.na(residuals(fit))), n_diff)
  expect_equal(sum(is.na(fitted(fit))), n_diff)
  expect_true(all(is.na(as.numeric(residuals(fit))[1:n_diff])))

  expect_equal(fit$nobs, length(y))
  expect_equal(fit$n.used, length(y) - n_diff)

  # fitted values are on the scale of the original series
  keep <- (n_diff + 1):length(y)
  expect_equal(
    as.numeric(fitted(fit))[keep] + as.numeric(residuals(fit))[keep],
    y[keep]
  )
})

test_that("no differencing leaves residuals unpadded", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  expect_equal(sum(is.na(residuals(fit))), 0)
  expect_equal(fit$n.used, fit$nobs)
  expect_equal(
    as.numeric(fitted(fit)) + as.numeric(residuals(fit)),
    y
  )
})

test_that("the start-up transient does not dominate the residuals", {
  y <- sim_sint(n = 300, s = 4)
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 1, 1)),
    seasonal.period = 4
  )
  r <- as.numeric(residuals(fit))
  r <- r[!is.na(r)]

  expect_lt(sd(r[1:20]) / sd(r), 2)
})

test_that("a series too short after differencing is rejected", {
  expect_error(
    msarima(
      sim_arma(n = 20),
      order = c(0, 0, 0),
      seasonal = list(c(0, 2, 0)),
      seasonal.period = 10
    ),
    "at least 10 are needed"
  )
})
