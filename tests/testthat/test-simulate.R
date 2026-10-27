test_that("msarima.sim returns a series of the requested length and class", {
  set.seed(3)
  y <- msarima.sim(
    nsim = 120,
    order = c(1, 0, 0),
    seasonal = list(c(1, 0, 0), c(0, 0, 1)),
    seasonal.period = c(4, 12),
    ar = c(0.5, 0.3),
    ma = 0.4,
    sigma2 = 1
  )

  expect_s3_class(y, "msts")
  expect_length(as.numeric(y), 120)
  expect_equal(attr(y, "msts"), c(4, 12))

  y1 <- msarima.sim(
    nsim = 50,
    order = c(1, 0, 0),
    seasonal = list(c(1, 0, 0)),
    seasonal.period = 4,
    ar = c(0.5, 0.3)
  )
  expect_s3_class(y1, "ts")
  expect_equal(frequency(y1), 4)
})

test_that("coefficient vectors must match the model orders", {
  expect_error(
    msarima.sim(
      nsim = 50,
      order = c(1, 0, 0),
      seasonal = list(c(1, 0, 0)),
      seasonal.period = 4,
      ar = 0.5
    ),
    "'ar' must have length 2"
  )
  expect_error(
    msarima.sim(
      nsim = 50,
      order = c(0, 0, 1),
      seasonal = list(c(0, 0, 1)),
      seasonal.period = 4,
      ma = 0.5
    ),
    "'ma' must have length 2"
  )
})

test_that("msarima.sim is centred on zero by default", {
  set.seed(5)
  y <- as.numeric(msarima.sim(
    nsim = 4000,
    order = c(1, 0, 0),
    seasonal = list(c(0, 0, 0)),
    seasonal.period = 4,
    ar = 0.5,
    sigma2 = 1
  ))

  expect_equal(mean(y), 0, tolerance = 0.15)
})

test_that("msarima.sim places the series at the level implied by intercept", {
  set.seed(5)
  y <- as.numeric(msarima.sim(
    nsim = 4000,
    order = c(1, 0, 0),
    seasonal = list(c(0, 0, 0)),
    seasonal.period = 4,
    ar = 0.5,
    sigma2 = 1,
    intercept = 25
  ))

  # the intercept is not the mean: the level is intercept / (1 - ar)
  expect_equal(mean(y), 50, tolerance = 0.5)
  # the level must be reached immediately, not built up during the series
  expect_equal(mean(y[1:200]), 50, tolerance = 2)
})

test_that("simulate.msarima reproduces the level of the fitted model", {
  set.seed(13)
  y <- 100 + as.numeric(msarima.sim(
    nsim = 300,
    order = c(1, 0, 0),
    seasonal = list(c(0, 0, 0)),
    seasonal.period = 4,
    ar = 0.5,
    sigma2 = 1
  ))
  fit <- msarima(
    y,
    order = c(1, 0, 0),
    seasonal = list(c(0, 0, 0)),
    seasonal.period = 4,
    include.constant = TRUE
  )

  sim <- as.numeric(simulate(fit, nsim = 2000, seed = 1))
  expect_equal(mean(sim), mean(y), tolerance = 0.05)
})

test_that("simulate.msarima is reproducible through 'seed'", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  a <- simulate(fit, nsim = 40, seed = 42)
  b <- simulate(fit, nsim = 40, seed = 42)

  expect_identical(as.numeric(a), as.numeric(b))
  expect_length(as.numeric(a), 40)
})
