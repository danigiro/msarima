test_that("print and summary work and return their input invisibly", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  expect_output(print(fit), "MSARIMA\\(1,0,1\\)\\(0,0,1\\)\\[4\\]")
  expect_output(print(fit), "Coefficients")
  capture.output(expect_invisible(print(fit)))

  s <- summary(fit)
  expect_s3_class(s, "summary_msarima")
  expect_equal(dim(s$coefficients), c(3L, 4L))
  expect_equal(
    colnames(s$coefficients),
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  )
  expect_equal(s$coefficients[, "Estimate"], coef(fit))
  expect_output(print(s), "Coefficients")
})

test_that("accessors return the stored components", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )

  expect_identical(coef(fit), fit$coef)
  expect_identical(fitted(fit), fit$fitted)
  expect_identical(residuals(fit), fit$residuals)
})

test_that("logLik, AIC and BIC are mutually consistent", {
  y <- sim_arma()
  fit <- msarima(
    y,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4
  )
  npar <- length(coef(fit)) + 1 # coefficients plus sigma

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), fit$loglik)
  expect_equal(attr(ll, "df"), npar)
  expect_equal(attr(ll, "nobs"), fit$n.used)

  expect_equal(AIC(fit), -2 * fit$loglik + 2 * npar)
  expect_equal(BIC(fit), -2 * fit$loglik + log(fit$n.used) * npar)
  expect_lt(AIC(fit), fit$aicc)
})
