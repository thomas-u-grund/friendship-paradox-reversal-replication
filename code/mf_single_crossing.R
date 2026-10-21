# Check of the single-crossing assumption in Proposition 3: Qcal > 1 on (0, lambda*), Qcal < 1 on (lambda*, nu*lambda*] (and beyond, up to large lambda).
mkS <- function(d, p) { p <- p / sum(p); mu1 <- sum(d * p); mu2 <- sum(d^2 * p); nu <- mu2/mu1 - 1
  H <- function(E) sapply(E, function(e) sum(p * (d - 1) * (-expm1(-d * e))) / mu1); list(d = d, p = p, mu1 = mu1, nu = nu, H = H) }
nbS <- function(mu, cv) { inv_r <- max(cv^2 - 1/mu, 1e-3); d <- 1:round(mu * (1 + 8 * cv)); mkS(d, dnbinom(d, size = 1/inv_r, mu = mu)) }
lnS <- function(mu, cv) { s2 <- log(1 + cv^2); m <- log(mu) - s2/2; d <- 1:round(mu * (1 + 12 * cv)); mkS(d, dlnorm(d, m, sqrt(s2))) }
plS <- function(mu, cv) { # truncated discrete power law with min degree dmin; exponent chosen to hit the mean, cutoff 30*mean
  dmax <- round(30 * mu); d <- 1:dmax; best <- NULL
  f <- function(a, dmin) { p <- ifelse(d >= dmin, d^(-a), 0); sum(d * p) / sum(p) }
  dmin <- max(2, round(mu / 3)); a <- tryCatch(uniroot(function(a) f(a, dmin) - mu, c(1.05, 6))$root, error = function(e) NA)
  if (is.na(a)) return(NULL); mkS(d, ifelse(d >= dmin, d^(-a), 0)) }
omega <- function(S, u, m = 60) { E <- u / S$nu^m; for (i in 1:m) E <- S$H(E); E }
Qcal <- function(S, lambda) sapply(lambda, function(lam) { num <- den <- 0
  for (i in seq_along(S$d)) { if (S$p[i] < 1e-6) next; s <- S$d[i]; Ek <- omega(S, lam * s / (S$mu1 * (S$nu - 1))); Ekm <- omega(S, (lam / S$nu) * s / (S$mu1 * (S$nu - 1)))
    pik <- exp(-S$d * Ekm) * (-expm1(-S$d * (Ek - Ekm))); num <- num + S$p[i] * sum(S$p * S$d * pik); den <- den + S$p[i] * sum(S$p * pik) }; num / (S$mu1 * den) })
check <- function(S, label) {
  lam <- 10^seq(-4, 3.5, length.out = 150); Q <- Qcal(S, lam); above <- Q >= 1
  nswitch <- sum(diff(above) != 0); j <- which(!above)[1]
  data.frame(law = label, nu = S$nu, ok_nu = S$nu > 1, Qsmall = Q[1], crossings = nswitch,
    single = (nswitch == 1) && above[1], lamstar = if (is.na(j)) NA else lam[j], minQ_after = if (is.na(j)) NA else min(Q[j:length(Q)]), maxQ_after = if (is.na(j)) NA else max(Q[j:length(Q)])) }
res <- list(); i <- 0
for (mu in c(6, 12, 25, 50)) for (cv in c(0.4, 1.0, 1.4)) {
  i <- i + 1; res[[i]] <- check(nbS(mu, cv), sprintf("NB mu=%g cv=%g", mu, cv))
  i <- i + 1; res[[i]] <- check(lnS(mu, cv), sprintf("LN mu=%g cv=%g", mu, cv)) }
for (mu in c(6, 12, 25)) { S <- plS(mu, 0); if (!is.null(S)) { i <- i + 1; res[[i]] <- check(S, sprintf("PL mu=%g", mu)) } }
R <- do.call(rbind, res); print(R, digits = 3, row.names = FALSE)
cat(sprintf("\n%d degree laws; nu>1 in %d; single crossing (Q>=1 then <1, exactly one switch) in %d; maximum of Q after the crossing: %.3f\n", nrow(R), sum(R$ok_nu), sum(R$single), max(R$maxQ_after, na.rm = TRUE)))
write.csv(R, "../results/mf_single_crossing.csv", row.names = FALSE)
