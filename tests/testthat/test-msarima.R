test_that("input is validated", {
  y <- sim_arma(n = 60)

  expect_error(msarima(letters), "must be numeric")
  expect_error(msarima(c(y[1:20], NA)), "missing values")
  expect_error(msarima(1:5), "at least 10 are needed")
  expect_error(
    msarima(y, seasonal = c(0, 0, 1), seasonal.period = 4),
    "must be a list"
  )
  expect_error(
    msarima(y, seasonal = list(c(0, 0, 1)), seasonal.period = c(4, 12)),
    "must have the same length"
  )
  expect_error(
    msarima(y, seasonal = list(c(0, 0)), seasonal.period = 4),
    "length 3"
  )
  expect_error(
    msarima(y, seasonal = list(c(0, 0, -1)), seasonal.period = 4),
    "non-negative"
  )
  # seasonal period must be >= 2 and at most half the series length
  expect_error(
    msarima(y, seasonal = list(c(1, 0, 0)), seasonal.period = 1),
    "Invalid 'seasonal.period"
  )
  expect_error(
    msarima(y, seasonal = list(c(1, 0, 0)), seasonal.period = 50),
    "Invalid 'seasonal.period"
  )
})

test_that("the fitted object has the documented structure", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  expect_s3_class(fit, "msarima")
  expect_equal(fit$method, "mSARIMA(1,0,1)(0,0,1)[4]")
  expect_equal(fit$order, c(1, 0, 1))
  expect_equal(fit$seasonal.period, 4)
  expect_equal(fit$series, "y")
  expect_true(is.matrix(fit$var_coef))
  expect_equal(dim(fit$var_coef), c(3L, 3L))
  expect_true(fit$sigma2 > 0)
})

test_that("coefficients are named and ordered AR, SAR, MA, SMA", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(2, 0, 1),
    seasonal = list(c(1, 0, 1), c(0, 0, 1)),
    seasonal.period = c(4, 12)
  )

  expect_named(
    coef(fit),
    c("AR(1)", "AR(2)", "SAR(1,1)", "MA(1)", "SMA(1,1)", "SMA(2,1)")
  )
  expect_equal(fit$method, "mSARIMA(2,0,1)(1,0,1)[4](0,0,1)[12]")
})

test_that("seasonal orders >= 2 get one coefficient per seasonal lag", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(0, 0, 0),
    seasonal = list(c(2, 0, 0)),
    seasonal.period = 4
  )

  expect_named(coef(fit), c("SAR(1,1)", "SAR(1,2)"))
})

test_that("the model can be fitted without any seasonal component", {
  y <- sim_arma()
  fit <- msarima(y, order = c(1, 0, 1), seasonal = list())

  expect_named(coef(fit), c("AR(1)", "MA(1)"))
  expect_equal(fit$method, "mSARIMA(1,0,1)")
})

test_that("known parameters are recovered on a long series", {
  set.seed(7)
  y <- as.numeric(msarima.sim(
    nsim = 800,
    order = c(1, 0, 0),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4,
    ar = 0.6,
    ma = -0.5,
    sigma2 = 1
  ))
  fit <- msarima(
    y,
    order = c(1, 0, 0),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  expect_equal(unname(coef(fit)), c(0.6, -0.5), tolerance = 0.1)
  expect_equal(fit$sigma2, 1, tolerance = 0.15)
})

test_that("fixed parameters are held at their value", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4,
    fixed = c(0.4, NA, NA, NA)
  )

  expect_equal(unname(coef(fit)[1]), 0.4)
  expect_equal(fit$mask, c(TRUE, FALSE, FALSE))

  # residuals must be the ones implied by the constrained parameters, not by
  # the unconstrained optimum
  free <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )
  expect_false(isTRUE(all.equal(
    as.numeric(residuals(fit)),
    as.numeric(residuals(free))
  )))
  expect_equal(
    as.numeric(fitted(fit)) + as.numeric(residuals(fit)),
    y
  )
})

test_that("'fixed' must have one entry per parameter plus sigma", {
  y <- sim_arma()
  expect_error(
    msarima(
      y,
      order = c(1, 0, 1),
      seasonal = list(c(0, 0, 1)),
      seasonal.period = 4,
      fixed = c(0.4, NA)
    ),
    "'fixed' has length 2"
  )
})

test_that("include.constant adds an intercept only without differencing", {
  y <- 20 + sim_arma()

  fit <- msarima(
    y,
    order = c(1, 0, 0),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4,
    include.constant = TRUE
  )
  expect_true("intercept" %in% names(coef(fit)))
  # the intercept is not the mean: the implied level is mu / (1 - phi)
  expect_equal(
    unname(coef(fit)["intercept"] / (1 - coef(fit)["AR(1)"])),
    mean(y),
    tolerance = 0.1
  )

  expect_warning(
    fit_d <- msarima(
      y,
      order = c(1, 1, 0),
      seasonal = list(c(0, 0, 1)),
      seasonal.period = 4,
      include.constant = TRUE
    ),
    "No constant term fitted"
  )
  expect_false("intercept" %in% names(coef(fit_d)))
})

test_that("ts and msts inputs supply the seasonal periods", {
  y <- sim_arma()

  fit_ts <- msarima(ts(y, frequency = 4), order = c(1, 0, 0), seasonal = list(c(1, 0, 0)))
  expect_equal(fit_ts$seasonal.period, 4)

  y_msts <- msts(y, seasonal.periods = c(4, 12))
  fit_msts <- msarima(
    y_msts,
    order = c(1, 0, 0),
    seasonal = list(c(1, 0, 0), c(0, 0, 1))
  )
  expect_equal(fit_msts$seasonal.period, c(4, 12))
})
