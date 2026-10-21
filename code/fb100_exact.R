# Exact/sampled Facebook100 analysis on the largest component of each network.
# Output: fb100_exact.csv  (one row per network). Run from any directory.
suppressMessages({library(igraph); library(parallel); library(dplyr); library(readr)})
root <- "../../"   # folder that contains facebook_networks/mtx_files/ (Facebook100 in Matrix Market format)
outdir <- "../results/"
files <- list.files(paste0(root,"facebook_networks/mtx_files/"), pattern="\\.mtx$", full.names=TRUE)
files <- files[file.info(files)$size < 15*1024*1024]
files <- files[order(-file.info(files)$size)]
K <- 6; NSRC <- 2000; EXACT_N <- 4000; NBOOT <- 200

load_lcc <- function(f) {
  l <- readLines(f); dl <- l[!grepl("^%", l)][-1]; e <- read.table(text=dl)[,1:2]; if (min(e)==0) e <- e+1
  g <- simplify(graph_from_edgelist(as.matrix(e), directed=FALSE))
  cc <- components(g); induced_subgraph(g, which(cc$membership == which.max(cc$csize)))
}
shellstats <- function(g, nsrc, K, seed) {
  set.seed(seed); n <- vcount(g); deg <- degree(g)
  src <- if (n <= EXACT_N) seq_len(n) else sample.int(n, nsrc)
  SZ <- matrix(0, length(src), K); DS <- SZ; ecc <- 0
  idx <- split(seq_along(src), ceiling(seq_along(src)/200))
  for (ii in idx) {
    D <- distances(g, v=src[ii]); ecc <- max(ecc, max(D))
    for (k in 1:K) { M <- (D == k); SZ[ii,k] <- rowSums(M); DS[ii,k] <- as.vector(M %*% deg) }
  }
  list(SZ=SZ, DS=DS, ecc=ecc, deg=deg, n=n, src=src, exact=(n<=EXACT_N))
}
summ <- function(o, mu1, mu2) {
  SZ <- o$SZ; DS <- o$DS; ns <- nrow(SZ)
  Bk <- function(idx) { b <- colSums(DS[idx,,drop=FALSE])/colSums(SZ[idx,,drop=FALSE]); b[1] <- mu2/mu1; b }
  B <- Bk(seq_len(ns)); S <- B/c(mu1, B[-length(B)])
  A <- sapply(1:ncol(SZ), function(k) { ok <- SZ[,k] > 0; mean(DS[ok,k]/SZ[ok,k]) })
  SA <- A/c(mu1, A[-length(A)]) ; SA[1] <- A[1]/mu1
  bs <- replicate(NBOOT, { i <- sample.int(ns, replace=TRUE); b <- Bk(i); b/c(mu1, b[-length(b)]) })
  list(B=B, S=S, A=A, SA=SA, se=apply(bs, 1, sd), En=colMeans(SZ))
}
mf_predict <- function(d, K=3) {
  n <- length(d); tab <- table(d); dv <- as.numeric(names(tab)); nd <- as.numeric(tab); S <- sum(d); U <- length(dv)
  size <- matrix(0,U,K); dsum <- matrix(0,U,K)
  for (s in seq_len(U)) { Fk <- dv[s]; cum <- 0; rprev <- rep(0,U)
    for (k in 1:K) { cum <- cum + Fk; r <- 1-exp(-dv*cum/S); newp <- pmax(r-rprev,0)
      size[s,k] <- sum(nd*newp); dsum[s,k] <- sum(nd*newp*dv); Fk <- sum(nd*newp*(dv-1)); rprev <- r } }
  colSums(nd*dsum)/colSums(nd*size)
}
nb_walk <- function(adj, deg, start, kmax) {
  cur <- start; prev <- 0L; out <- rep(NA_real_, kmax)
  for (s in seq_len(kmax)) { nb <- adj[[cur]]; if (prev != 0L) nb <- nb[nb != prev]; if (!length(nb)) break
    nxt <- if (length(nb)==1) nb else nb[sample.int(length(nb),1)]; prev <- cur; cur <- nxt; out[s] <- deg[cur] }
  out
}
one <- function(f) {
  t0 <- Sys.time(); nm <- basename(f); g <- load_lcc(f); n <- vcount(g); m <- ecount(g); deg <- degree(g)
  mu1 <- mean(deg); mu2 <- mean(deg^2)
  o <- shellstats(g, NSRC, K, seed=1); r <- summ(o, mu1, mu2)
  # walk statistic (non-backtracking), start vertices sampled, 10 walks each
  set.seed(2); adj <- as_adj_list(g, mode="all"); st <- if (n<=EXACT_N) seq_len(n) else sample.int(n, NSRC)
  wk <- matrix(NA_real_, length(st)*10, 3); row <- 1L
  for (rep in 1:10) for (s in st) { wk[row,] <- nb_walk(adj, deg, s, 3); row <- row+1L }
  Wk <- colMeans(wk, na.rm=TRUE); Sw <- Wk/c(mu1, Wk[-3])
  mf <- mf_predict(deg, 3)
  # degree-preserving null
  set.seed(3); g0 <- rewire(g, keeping_degseq(niter=10*m)); cc0 <- components(g0)
  g0 <- induced_subgraph(g0, which(cc0$membership==which.max(cc0$csize)))
  d0 <- degree(g0); o0 <- shellstats(g0, 1000, 3, seed=4); r0 <- summ(o0, mean(d0), mean(d0^2))
  data.frame(network=nm, n=n, m=m, mu1=mu1, cv=sd(deg)/mu1, maxdeg=max(deg),
    clustering=transitivity(g, type="global"), assort=assortativity_degree(g), ecc_max=o$ecc, exact=o$exact,
    t(setNames(r$S[1:K], paste0("S",1:K))), t(setNames(r$se[1:K], paste0("seS",1:K))),
    t(setNames(r$B/mu1, paste0("Q",1:K))), t(setNames(r$A/mu1, paste0("Aq",1:K))),
    t(setNames(r$SA[1:3], paste0("SA",1:3))), t(setNames(r$En/n, paste0("cov",1:K))),
    t(setNames(Sw, paste0("Sw",1:3))), t(setNames(mf/mu1, paste0("MFq",1:3))),
    null_S2=r0$S[2], null_S3=r0$S[3], null_Q3=r0$B[3]/mean(d0), secs=as.numeric(Sys.time()-t0, units="secs"))
}
res <- mclapply(files, function(f) tryCatch({ x <- one(f); cat(basename(f), round(x$secs), "s\n", file=paste0(outdir,"fb100_exact.log"), append=TRUE); x },
         error=function(e) { cat("FAILED", basename(f), conditionMessage(e), "\n", file=paste0(outdir,"fb100_exact.log"), append=TRUE); NULL }),
         mc.cores=min(4, detectCores()-1), mc.preschedule=FALSE)
write_csv(bind_rows(res), paste0(outdir,"fb100_exact.csv"))
cat("DONE\n", file=paste0(outdir,"fb100_exact.log"), append=TRUE)
