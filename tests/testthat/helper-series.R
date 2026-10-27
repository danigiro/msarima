# Small simulated series shared across test files

# Stationary SARMA(1,0,1)(0,0,1)[4]
sim_arma <- function(n = 150, seed = 101) {
  set.seed(seed)
  as.numeric(msarima.sim(
    nsim = n,
    order = c(1, 0, 1),
    seasonal = list(c(0, 0, 1)),
    seasonal.period = 4,
    ar = 0.5,
    ma = c(0.3, -0.4),
    sigma2 = 1
  ))
}

# Same process integrated once, i.e. an ARIMA(1,1,1)(0,0,1)[4] around level 50
sim_int <- function(n = 150, seed = 101) {
  50 + cumsum(sim_arma(n = n, seed = seed))
}

# Same process with a seasonal unit root, i.e. (0,1,1)[4]
sim_sint <- function(n = 160, seed = 202, s = 4) {
  z <- sim_arma(n = n, seed = seed)
  y <- numeric(n)
  y[1:s] <- z[1:s]
  for (i in (s + 1):n) {
    y[i] <- y[i - s] + z[i]
  }
  y
}
