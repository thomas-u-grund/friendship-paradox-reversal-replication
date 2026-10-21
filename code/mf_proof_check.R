# Check of the 1D reduction E_{k+1} = delta_s + H(E_k) and of the Koenigs/Poincare form Q_k(n) = Qcal(nu^k / n).
nbpmf <- function(mu, cv) { inv_r <- max(cv^2 - 1/mu, 1e-3); d <- 1:round(mu * (1 + 8 * cv)); p <- dnbinom(d, size = 1/inv_r, mu = mu); list(d = d, p = p / sum(p)) }
setup <- function(mu, cv) { z <- nbpmf(mu, cv); d <- z$d; p <- z$p; mu1 <- sum(d * p); mu2 <- sum(d^2 * p); nu <- mu2/mu1 - 1
  H <- function(E) sapply(E, function(e) sum(p * (d - 1) * (-expm1(-d * e))) / mu1); list(d = d, p = p, mu1 = mu1, nu = nu, H = H) }
# direct recursion in (E, x) as in the paper vs the 1D reduction
direct_E <- function(S, s, n, K) { E <- 0; x <- s / (n * S$mu1); out <- numeric(K); for (k in 1:K) { Enew <- E + x; pi_k <- exp(-S$d * E) * (-expm1(-S$d * (Enew - E))); x <- sum(S$p * (S$d - 1) * pi_k) / S$mu1; E <- Enew; out[k] <- E }; out }
oned_E <- function(S, s, n, K) { delta <- s / (n * S$mu1); E <- 0; out <- numeric(K); for (k in 1:K) { E <- delta + S$H(E); out[k] <- E }; out }
S <- setup(12, 1.0); a <- direct_E(S, 12, 1e6, 12); b <- oned_E(S, 12, 1e6, 12); cat("max |E_k(direct) - E_k(1D reduction)| over k=1..12:", format(max(abs(a - b)), digits = 3), "\n")
omega <- function(S, u, m = 60) { E <- u / S$nu^m; for (i in 1:m) E <- S$H(E); E }          # Poincare function: H(omega(u)) = omega(nu u)
chk <- max(abs(S$H(omega(S, c(0.01, 0.1, 1, 10))) - omega(S, S$nu * c(0.01, 0.1, 1, 10)))); cat("max |H(omega(u)) - omega(nu u)|:", format(chk, digits = 3), "\n")
Qcal <- function(S, lambda) sapply(lambda, function(lam) { num <- den <- 0
  for (i in seq_along(S$d)) { s <- S$d[i]; Ek <- omega(S, lam * s / (S$mu1 * (S$nu - 1))); Ekm <- omega(S, (lam / S$nu) * s / (S$mu1 * (S$nu - 1))); pik <- exp(-S$d * Ekm) * (-expm1(-S$d * (Ek - Ekm)))
    num <- num + S$p[i] * sum(S$p * S$d * pik); den <- den + S$p[i] * sum(S$p * pik) }; num / (S$mu1 * den) })
# compare Q_k(n) from the recursion with Qcal(nu^k/n)
Qk_direct <- function(S, n, K = 30) { num <- den <- numeric(K); for (i in seq_along(S$d)) { E <- oned_E(S, S$d[i], n, K); Em <- c(0, E[-K]); pik <- sapply(1:K, function(k) exp(-S$d * Em[k]) * (-expm1(-S$d * (E[k] - Em[k])))); num <- num + S$p[i] * colSums(S$p * S$d * pik); den <- den + S$p[i] * colSums(S$p * pik) }; num / (S$mu1 * den) }
cat('Qcal at lambda = nu^k/n (n=1e5), k=1..8:', round(Qcal(S, S$nu^(1:8)/1e5), 4), '\n'); cat('direct Q_k(n=1e5), k=1..8:', round(Qk_direct(S, 1e5)[1:8], 4), '\n')
lam <- 10^seq(-3, 4, length.out = 400); Q <- Qcal(S, lam); j <- which(Q < 1)[1]; lamstar <- exp(approx(Q[(j-1):j], log(lam[(j-1):j]), xout = 1)$y)
cat(range(Q), '\n'); cat(sprintf("NB 12/1.0: nu = %.3f; Qcal(lambda) at small lambda = %.3f (mu2/mu1^2 = %.3f); lambda* = %.4f; c* = log_nu(lambda*) = %.3f\n", S$nu, Q[1], sum(S$d^2*S$p)/S$mu1^2, lamstar, log(lamstar)/log(S$nu)))
# intrinsic crossing k*(n) = log_nu(n lambda*) vs interpolated crossing from the recursion
cross <- function(Q) { j <- which(Q[-1] < 1)[1]; k <- j + 1; (k - 1) + (Q[k-1] - 1) / (Q[k-1] - Q[k]) }
ns <- 10^seq(3, 9, by = 0.25); kt <- sapply(ns, function(n) cross(Qk_direct(S, n))); kstar <- log(ns * lamstar) / log(S$nu); d <- kt - kstar
cat(sprintf("interpolated crossing minus intrinsic crossing over n=1e3..1e9: mean %.3f, range %.3f to %.3f (periodic in log_nu n)\n", mean(d), min(d), max(d)))

for (n in c(1e3, 1e5, 1e7)) { Qd <- Qk_direct(S, n); k <- 2:8; Qc <- Qcal(S, S$nu^k / n); cat(sprintf("n=1e%d: max_k |Q_k(n) - Qcal(nu^k/n)| over k=2..8: %.2e\n", log10(n), max(abs(Qd[k] - Qc)))) }
