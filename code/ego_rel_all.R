# Individual-level paradox at distance k: share of egos whose distance-k shell has higher mean degree than the ego.
src <- readLines("fb100_exact.R"); i <- grep("^res <- mclapply", src); eval(parse(text=src[1:(i-1)]))
K5 <- 5
mf_class <- function(d, K) {              # expected shell mean degree for each source degree class
  n <- length(d); tab <- table(d); dv <- as.numeric(names(tab)); nd <- as.numeric(tab); S <- sum(d); U <- length(dv)
  m <- matrix(NA_real_, U, K)
  for (s in seq_len(U)) { Fk <- dv[s]; cum <- 0; rprev <- rep(0,U)
    for (k in 1:K) { cum <- cum + Fk; r <- 1-exp(-dv*cum/S); newp <- pmax(r-rprev,0)
      sz <- sum(nd*newp); if (sz > 1e-6) m[s,k] <- sum(nd*newp*dv)/sz
      Fk <- sum(nd*newp*(dv-1)); rprev <- r } }
  list(dv=dv, nd=nd, m=m)
}
one_ego <- function(f) {
  g <- load_lcc(f); n <- vcount(g); deg <- degree(g); mu1 <- mean(deg)
  set.seed(7); egos <- if (n <= EXACT_N) seq_len(n) else sample.int(n, NSRC)
  M <- matrix(NA_real_, length(egos), K5)
  idx <- split(seq_along(egos), ceiling(seq_along(egos)/200))
  for (ii in idx) { D <- distances(g, v=egos[ii])
    for (k in 1:K5) { Mk <- (D==k); sz <- rowSums(Mk); ds <- as.vector(Mk %*% deg); M[ii,k] <- ifelse(sz>0, ds/sz, NA) } }
  de <- deg[egos]; out <- list(network=basename(f), n=n, mu1=mu1)
  terc <- cut(de, quantile(de, c(0,1/3,2/3,1)), include.lowest=TRUE, labels=FALSE)
  for (k in 1:K5) { ok <- !is.na(M[,k]); out[[paste0("own",k)]] <- mean(M[ok,k] > de[ok]); out[[paste0("pop",k)]] <- mean(M[ok,k] > mu1)
    out[[paste0("ratio",k)]] <- mean(M[ok,k]/de[ok]); out[[paste0("covered",k)]] <- mean(ok)
    for (tt in 1:3) out[[paste0("own",k,"_t",tt)]] <- if (k<=4) mean(M[ok & terc==tt,k] > de[ok & terc==tt]) else NA }
  mc <- mf_class(deg, 4)
  for (k in 1:4) { ok <- !is.na(mc$m[,k]); out[[paste0("mf_own",k)]] <- sum(mc$nd[ok & mc$m[,k] > mc$dv])/sum(mc$nd[ok])
    out[[paste0("mf_pop",k)]] <- sum(mc$nd[ok & mc$m[,k] > mu1])/sum(mc$nd[ok]) }
  as.data.frame(out)
}
res <- mclapply(files, function(f) tryCatch(one_ego(f), error=function(e) {cat("FAILED", basename(f), conditionMessage(e), "\n"); NULL}),
                mc.cores=min(4, detectCores()-1), mc.preschedule=FALSE)
write_csv(bind_rows(res), "ego_rel_all.csv"); cat("DONE\n")
