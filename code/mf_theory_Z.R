# WP1 theory pilots on the deterministic mean-field model (continuous class counts n_d = n p_d).
suppressMessages(library(parallel))
mf_pmf <- function(dv, pd, n, K = 60) {
  nd <- n * pd; S <- sum(nd * dv); mu1 <- S / n; U <- length(dv); size <- dsum <- matrix(0, U, K)
  for (s in seq_len(U)) { Fk <- dv[s]; cum <- 0; rprev <- rep(0, U)
    for (k in 1:K) { cum <- cum + Fk; r <- 1 - exp(-dv * cum / S); newp <- pmax(r - rprev, 0); size[s, k] <- sum(nd * newp); dsum[s, k] <- sum(nd * newp * dv); Fk <- sum(nd * newp * (dv - 1)); rprev <- r; if (Fk < 1e-12) break } }
  m <- dsum / pmax(size, 1e-300); m[size < 1e-9] <- NA
  P <- sapply(1:K, function(k) { ok <- !is.na(m[, k]); if (!any(ok)) NA else sum(nd[ok & m[, k] > dv]) / sum(nd[ok]) })
  list(Q = colSums(nd * dsum) / pmax(colSums(nd * size), 1e-300) / mu1, P = P) }
cross <- function(Q) { j <- which(Q[-1] < 1)[1]; if (is.na(j)) NA else { k <- j + 1; (k - 1) + (Q[k - 1] - 1) / (Q[k - 1] - Q[k]) } }
k50 <- function(P) { j <- which(P < 0.5)[1]; if (is.na(j)) NA else if (j == 1) 1 else (j - 1) + (P[j - 1] - 0.5) / (P[j - 1] - P[j]) }
pmf_family <- function(fam, mu, cv, dmax) { d <- 1:dmax
  p <- switch(fam,
    NB = { inv_r <- max(cv^2 - 1/mu, 1e-3); dnbinom(d, size = 1/inv_r, mu = mu) },
    LN = { s2 <- log(1 + cv^2); plnorm(d + 0.5, log(mu) - s2/2, sqrt(s2)) - plnorm(d - 0.5, log(mu) - s2/2, sqrt(s2)) },
    PL = { best <- NULL; for (dmin in 1:floor(mu)) { dd <- dmin:dmax; pp <- dd^-2.5; pp <- pp/sum(pp); m <- sum(dd*pp); if (is.null(best) || abs(m - mu) < best$e) best <- list(e = abs(m - mu), dmin = dmin) }
           ifelse(d >= best$dmin, d^-2.5, 0) })
  p <- p / sum(p); list(d = d, p = p) }
ns <- 10^seq(3, 12, by = 1)
# E1: scaling in other families / larger n
e1 <- expand.grid(fam = c("NB", "LN", "PL"), mu = c(12, 50), cv = c(0.7, 1.4), stringsAsFactors = FALSE); e1 <- e1[!(e1$fam == "PL" & e1$cv == 0.7), ]; e1$cv[e1$fam == "PL"] <- NA
out1 <- mclapply(seq_len(nrow(e1)), function(i) { f <- e1[i, ]; dmax <- if (f$fam == "PL") 2000 else round(f$mu * (1 + 6 * ifelse(is.na(f$cv), 1, f$cv)))
  dd <- pmf_family(f$fam, f$mu, f$cv, dmax); mu1 <- sum(dd$d * dd$p); mu2 <- sum(dd$d^2 * dd$p); nu <- mu2 / mu1 - 1
  kt <- sapply(ns, function(n) cross(mf_pmf(dd$d, dd$p, n)$Q)); z <- log(ns) / log(nu); fit <- lm(kt ~ z)
  sprintf("E1 %s mu=%d cv=%s: realised mean %.1f cv %.2f nu %.1f | slope %.3f intercept %.3f resid sd %.3f | crossing n=1e3 %.2f n=1e12 %.2f", f$fam, f$mu, ifelse(is.na(f$cv), "-", f$cv), mu1, sqrt(mu2 - mu1^2) / mu1, nu, coef(fit)[2], coef(fit)[1], sd(resid(fit)), kt[1], kt[length(kt)]) }, mc.cores = 4)
cat(unlist(out1), sep = "\n")
# E2: sign of the heterogeneity effect at fixed n = 1e4, by mean degree
e2 <- expand.grid(mu = c(6, 12, 25, 50, 100, 200), cv = c(0.4, 0.7, 1.0, 1.4, 2.0))
out2 <- mclapply(seq_len(nrow(e2)), function(i) { f <- e2[i, ]; dd <- pmf_family("NB", f$mu, f$cv, round(f$mu * (1 + 8 * f$cv))); kt <- cross(mf_pmf(dd$d, dd$p, 1e4)$Q); c(mu = f$mu, cv = f$cv, ktil = kt) }, mc.cores = 4)
r2 <- do.call(rbind, out2); cat("\nE2: continuous crossing at n = 1e4 (rows mean degree, columns CV 0.4,0.7,1.0,1.4,2.0)\n")
for (m in sort(unique(r2[, "mu"]))) { v <- r2[r2[, "mu"] == m, "ktil"]; cat(sprintf("mu=%3d: %s   change over CV 0.4->2.0: %+.2f\n", m, paste(sprintf("%.2f", v), collapse = " "), v[5] - v[1])) }
# E3: individual-level 50% crossing vs reversal distance, NB cv=1, mu=12/50
for (mu in c(12, 50)) { dd <- pmf_family("NB", mu, 1.0, round(mu * 9)); res <- sapply(c(1e3, 1e5, 1e7, 1e9), function(n) { o <- mf_pmf(dd$d, dd$p, n); c(ktil = cross(o$Q), kP50 = k50(o$P)) }); cat(sprintf("\nE3 NB cv=1 mu=%d: n=1e3,1e5,1e7,1e9\n  ktil: %s\n  kP50: %s\n", mu, paste(sprintf("%.2f", res["ktil", ]), collapse = " "), paste(sprintf("%.2f", res["kP50", ]), collapse = " "))) }
