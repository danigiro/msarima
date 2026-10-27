test_that("seasonal period order does not matter", {
  set.seed(1)
  z <- as.numeric(msarima.sim(600, seasonal = list(c(1,0,0), c(2,0,0)),
                              seasonal.period = c(4, 12), ar = c(0.4, 0.5, 0.3)))
  f1 <- msarima(z, seasonal = list(c(1,0,0), c(2,0,0)), seasonal.period = c(4, 12))
  f2 <- msarima(z, seasonal = list(c(2,0,0), c(1,0,0)), seasonal.period = c(12, 4))
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-6)
  expect_equal(unname(coef(f1)), unname(coef(f2)[c(3, 1, 2)]), tolerance = 1e-4)
  expect_equal(as.numeric(forecast(f1, h = 24)$mean),
               as.numeric(forecast(f2, h = 24)$mean), tolerance = 1e-4)
})

test_that("seasonal period order does not matter in simulation", {
  set.seed(2)
  s1 <- msarima.sim(300, seasonal = list(c(1,0,0), c(2,0,0)),
                    seasonal.period = c(4, 12), ar = c(.4, .5, .3))
  set.seed(2)
  s2 <- msarima.sim(300, seasonal = list(c(2,0,0), c(1,0,0)),
                    seasonal.period = c(12, 4), ar = c(.5, .3, .4))
  expect_equal(as.numeric(s1), as.numeric(s2))
})

test_that("duplicated seasonal periods are rejected", {
  y <- rnorm(200)
  expect_error(msarima(y, seasonal = list(c(1,0,0), c(1,0,0)), seasonal.period = c(12, 12)), "duplicated")
  expect_error(msarima.sim(200, seasonal = list(c(1,0,0), c(1,0,0)), seasonal.period = c(12, 12), ar = c(.3,.3)), "duplicated")
})
