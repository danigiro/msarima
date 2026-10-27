# Reference implementation of the original two-season Ljung-Box test, used
# below to check msportmanteau() reproduces it exactly for the same lag set.
ref_ljung_box_2s <- function(x, s1, s2, m1, m2, n.para) {
  n <- length(x)
  L <- sort(unique(c(s1 * (1:m1), s2 * (1:m2))))
  obs <- as.vector(acf(x, max(L), plot = FALSE)$acf^2)[L + 1]
  STATISTIC <- n * (n + 2) * sum(obs / (n - L))
  PARAMETER <- length(L) - n.para
  list(statistic = STATISTIC, df = PARAMETER,
       pvalue = pchisq(STATISTIC, PARAMETER, lower.tail = FALSE))
}

test_that("msportmanteau reproduces the original hand-rolled two-season formula", {
  set.seed(42)
  x <- rnorm(300)
  ref <- ref_ljung_box_2s(x, s1 = 4, s2 = 11, m1 = 5, m2 = 3, n.para = 1)

  out <- msportmanteau(x, season = c(4, 11), nlag = c(5, 3), npar = 1)

  expect_equal(unname(out$statistic), ref$statistic)
  expect_equal(unname(out$parameter), ref$df)
  expect_equal(out$p.value, ref$pvalue)
  expect_equal(out$lag, sort(unique(c(4 * (1:5), 11 * (1:3)))))
})

test_that("lags_short adds the low-order non-seasonal lags", {
  x <- rnorm(200)
  out <- msportmanteau(x, season = 4, nlag = 2, lags_short = 3)
  expect_equal(out$lag, sort(unique(c(1:3, 4 * (1:2)))))
})

test_that("white noise gives a large p-value, a seasonal series a small one", {
  set.seed(1)
  wn <- rnorm(400)
  expect_gt(msportmanteau(wn, season = 12, nlag = 3)$p.value, 0.05)
  expect_gt(msportmanteau(wn, type = "Monti", season = 12, nlag = 3)$p.value, 0.05)

  s <- 12
  seasonal <- rep(rnorm(s) * 5, length.out = 400) + rnorm(400, sd = 0.1)
  expect_lt(msportmanteau(seasonal, season = s, nlag = 3)$p.value, 0.001)
  expect_lt(msportmanteau(seasonal, type = "Monti", season = s, nlag = 3)$p.value, 0.001)
})

test_that("npar reduces the degrees of freedom", {
  x <- rnorm(200)
  out0 <- msportmanteau(x, season = 4, nlag = 3)
  out2 <- msportmanteau(x, season = 4, nlag = 3, npar = 2)
  expect_equal(unname(out0$parameter), length(out0$lag))
  expect_equal(unname(out2$parameter), length(out2$lag) - 2)
})

test_that("type = \"Monti\" uses the partial autocorrelations, not the acf", {
  set.seed(5)
  x <- rnorm(200)
  b <- msportmanteau(x, season = 4, nlag = 3)
  m <- msportmanteau(x, type = "Monti", season = 4, nlag = 3)
  expect_false(isTRUE(all.equal(unname(b$statistic), unname(m$statistic))))
  expect_equal(b$lag, m$lag)
})

test_that("season and npar are picked up from a fitted msarima object", {
  y <- sim_arma()
  fit <- msarima(y, order = c(1, 0, 1), seasonal = list(c(0, 0, 1)),
                 seasonal.period = 4)
  # model without differencing; a differenced one is tested below
  out <- msportmanteau(fit, nlag = 5)
  expect_equal(out$lag, sort(unique(4 * (1:5))))
  expect_equal(unname(out$parameter), length(out$lag) - sum(!fit$mask))
  expect_equal(out$data.name, fit$series)
})

test_that("msarima residuals with leading NAs (differencing) are handled", {
  y <- sim_sint(n = 200, s = 4)
  fit <- msarima(y, order = c(0, 0, 0), seasonal = list(c(0, 1, 1)),
                 seasonal.period = 4)
  expect_true(anyNA(fit$residuals))

  out <- msportmanteau(fit, nlag = 3)
  n_clean <- length(na.contiguous(as.numeric(fit$residuals)))
  expect_lt(max(out$lag), n_clean)
  expect_false(is.na(out$statistic))
  expect_false(is.na(out$p.value))
})

test_that("empty season errors", {
  expect_error(msportmanteau(rnorm(50), season = numeric(0)), "season")
  expect_error(msportmanteau(rnorm(50), season = c(NA, 1)), "season")
})

test_that("non-positive degrees of freedom errors", {
  expect_error(msportmanteau(rnorm(50), season = 4, nlag = 1, npar = 10), "npar")
})

test_that("lags longer than the usable series error", {
  expect_error(msportmanteau(rnorm(20), season = 4, nlag = 10), "observations")
})

test_that("season is taken from an msts object when not supplied", {
  y <- msts(rnorm(200), seasonal.periods = c(4, 11))
  out <- msportmanteau(y, nlag = 2)
  expect_equal(out$lag, sort(unique(c(4 * (1:2), 11 * (1:2)))))
})

test_that("result has class htest with the expected fields", {
  out <- msportmanteau(rnorm(100), season = 4, nlag = 2)
  expect_s3_class(out, "htest")
  expect_named(out$statistic, "X-squared")
  expect_named(out$parameter, "df")
  expect_type(out$method, "character")
})
