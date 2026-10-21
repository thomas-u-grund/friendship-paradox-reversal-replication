# Mean-field check of the threshold structure for k>=3: is m_k(s) - s strictly decreasing in start degree s (single crossing of zero)? Is m_k(s) non-increasing?
src <- readLines("mf_single_crossing.R"); src <- src[1:grep("^omega <- ", src)-0]   # helper definitions only (mkS, nbS, lnS, plS, omega)
eval(parse(text = src[1:(grep("^omega <- ", src))]))
oned <- function(S, s, n, K) { delta <- s / (n * S$mu1); E <- 0; out <- numeric(K); for (k in 1:K) { E <- delta + S$H(E); out[k] <- E }; out }
shellmean <- function(S, n, K = 6) { keep <- which(S$p > 1e-8); M <- matrix(NA, length(keep), K)
  for (j in seq_along(keep)) { E <- oned(S, S$d[keep[j]], n, K); Em <- c(0, E[-K])
    for (k in 1:K) { pik <- exp(-S$d * Em[k]) * (-expm1(-S$d * (E[k] - Em[k]))); if (sum(S$p * pik) > 1e-10) M[j, k] <- sum(S$p * S$d * pik) / sum(S$p * pik) } }
  list(s = S$d[keep], M = M) }
laws <- list(); for (mu in c(6, 12, 25, 50)) for (cv in c(0.4, 1.0, 1.4)) { laws[[length(laws) + 1]] <- list(nbS(mu, cv), sprintf("NB %g/%g", mu, cv)); laws[[length(laws) + 1]] <- list(lnS(mu, cv), sprintf("LN %g/%g", mu, cv)) }
for (mu in c(6, 12, 25)) { S <- plS(mu, 0); if (!is.null(S)) laws[[length(laws) + 1]] <- list(S, sprintf("PL %g", mu)) }
res <- NULL
for (L in laws) for (n in c(1e3, 1e4, 1e5)) { r <- shellmean(L[[1]], n); for (k in 1:5) { m <- r$M[, k]; ok <- !is.na(m); if (sum(ok) < 3) next; d1 <- diff(m[ok]); d2 <- diff(m[ok] - r$s[ok])
  res <- rbind(res, data.frame(law = L[[2]], n = n, k = k, nonincr = all(d1 <= 1e-9), strict_dec_gap = all(d2 < 0), maxup = max(d1))) } }
print(aggregate(cbind(nonincr, strict_dec_gap) ~ k, res, function(x) paste0(sum(x), "/", length(x))))
print(res[!res$nonincr, ], digits = 3); write.csv(res, "../results/mf_threshold_check.csv", row.names = FALSE)
