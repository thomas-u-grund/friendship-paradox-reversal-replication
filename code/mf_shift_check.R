# Test of the translation property of the mean-field recursion: Q_{k+1}(nu * n) ~ Q_k(n) as n grows.
# Autonomous form: with E_k = cumulative exposure/Sigma and x_k = F_k/Sigma,  x_{k+1} = (1/mu1) sum_d p_d (d-1) (exp(-d E_{k-1}) - exp(-d E_k)), E_k = E_{k-1} + x_k, started at (E_0, x_1) = (0, s/(n mu1)).
mf_Q <- function(dv, pd, n, K = 80) {
  mu1 <- sum(dv * pd); U <- length(dv); num <- den <- numeric(K)
  for (s in seq_len(U)) { x <- dv[s] / (n * mu1); E <- 0
    for (k in 1:K) { Enew <- E + x; pi_k <- exp(-dv * E) - exp(-dv * Enew); num[k] <- num[k] + pd[s] * sum(pd * dv * pi_k); den[k] <- den[k] + pd[s] * sum(pd * pi_k)
      x <- sum(pd * (dv - 1) * pi_k) / mu1; E <- Enew; if (x < 1e-300) break } }
  list(Q = num / pmax(den, 1e-300) / mu1, share = den) }   # share[k] = expected shell-k fraction of the network
for (cfg in list(c(mu = 12, cv = 1.0), c(mu = 50, cv = 1.0), c(mu = 12, cv = 1.4))) {
  mu <- cfg["mu"]; cv <- cfg["cv"]; inv_r <- max(cv^2 - 1/mu, 1e-3); dv <- 1:round(mu * (1 + 8 * cv)); pd <- dnbinom(dv, size = 1/inv_r, mu = mu); pd <- pd / sum(pd)
  mu1 <- sum(dv * pd); nu <- sum(dv^2 * pd) / mu1 - 1; cat(sprintf("\nNB mean degree %d CV %.1f: nu = %.2f\n", mu, cv, nu))
  for (n0 in 10^c(3, 5, 7, 9, 11)) { a <- mf_Q(dv, pd, n0); b <- mf_Q(dv, pd, nu * n0); k <- 2:60; ok <- a$share[k] > 1e-8 & a$share[k] < 0.999 & !is.na(a$Q[k]) & !is.na(b$Q[k + 1])
    d <- abs(b$Q[k + 1] - a$Q[k])[ok]; cat(sprintf("  n = 1e%d: max |Q_{k+1}(nu n) - Q_k(n)| over shells holding 1e-8..0.999 of the network: %.2e  (%d shells)\n", log10(n0), max(d), length(d))) } }
