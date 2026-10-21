suppressMessages(library(parallel))
mf_profile <- function(dv, pd, n, K=40) {
  nd <- n*pd; S <- sum(nd*dv); mu1 <- S/n; U <- length(dv); size <- dsum <- matrix(0, U, K)
  for (s in seq_len(U)) { Fk <- dv[s]; cum <- 0; rprev <- rep(0,U)
    for (k in 1:K) { cum <- cum + Fk; r <- 1-exp(-dv*cum/S); newp <- pmax(r-rprev,0); size[s,k] <- sum(nd*newp); dsum[s,k] <- sum(nd*newp*dv); Fk <- sum(nd*newp*(dv-1)); rprev <- r; if (Fk < 1e-12) break } }
  colSums(nd*dsum)/pmax(colSums(nd*size),1e-300)/mu1 }
cross <- function(Q) { k <- which(Q[2:length(Q)] < 1)[1]+1; if (is.na(k)) NA else k-1+(Q[k-1]-1)/(Q[k-1]-Q[k]) }
grid <- expand.grid(cv=c(0.6,1.0,1.4), mu=c(12,50))
out <- mclapply(seq_len(nrow(grid)), function(i) { cvT <- grid$cv[i]; mu <- grid$mu[i]
  inv_r <- max(cvT^2-1/mu,1e-3); dmax <- round(mu*(1+6*cvT)); d <- 1:dmax; p <- dnbinom(d, size=1/inv_r, mu=mu); p <- p/sum(p)
  mu1 <- sum(d*p); mu2 <- sum(d^2*p); nu <- mu2/mu1-1; ns <- 10^seq(3,9,by=1); kt <- sapply(ns, function(n) cross(mf_profile(d,p,n)))
  z <- log(ns)/log(nu); f <- lm(kt~z)
  sprintf("CV %.1f mean deg %d: nu=%.1f  slope=%.3f  intercept=%.3f  resid sd=%.3f | k_cross at n=1e3..1e9: %s", cvT, mu, nu, coef(f)[2], coef(f)[1], sd(resid(f)), paste(round(kt,2), collapse=" ")) }, mc.cores=6)
cat(unlist(out), sep="\n")
