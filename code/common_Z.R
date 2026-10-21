# Version Z helpers (WP0). Definitions follow ANALYSIS_PLAN.md.
#   Q_k = B_k / mu1,  B_k = sum_v d_v n_k(v) / sum_v n_k(v)  (pooled enrichment, Theorem 1(a); B_1 = mu2/mu1 exactly)
#   k_rev  = min{k >= 2 : Q_k < 1}                         (discrete reversal distance)
#   ktil   = linear interpolation of the crossing of Q = 1  (continuous reversal distance, primary outcome)
#   P_k    = share of egos with n_k > 0 whose shell mean degree exceeds their own degree
suppressMessages({library(igraph); library(parallel); library(dplyr); library(readr)})

lcc <- function(g) { cc <- components(g); induced_subgraph(g, which(cc$membership == which.max(cc$csize))) }

# Shell profile from a set of source vertices. Returns per-source shell sizes SZ and shell degree sums DS (sources x K), K = largest distance seen.
shell_profile <- function(g, egos = NULL, chunk = 200) {
  n <- vcount(g); deg <- degree(g); if (is.null(egos)) egos <- seq_len(n)
  idx <- split(seq_along(egos), ceiling(seq_along(egos) / chunk))
  parts <- lapply(idx, function(ii) {
    D <- distances(g, v = egos[ii]); kmax <- max(D[is.finite(D)])
    sz <- ds <- matrix(0, length(ii), kmax)
    for (k in seq_len(kmax)) { Mk <- (D == k); sz[, k] <- rowSums(Mk); ds[, k] <- as.vector(Mk %*% deg) }
    list(sz = sz, ds = ds, kmax = kmax) })
  K <- max(sapply(parts, `[[`, "kmax")); pad <- function(m) cbind(m, matrix(0, nrow(m), K - ncol(m)))
  list(SZ = do.call(rbind, lapply(parts, function(p) pad(p$sz))), DS = do.call(rbind, lapply(parts, function(p) pad(p$ds))),
       deg = deg, n = n, egos = egos, exact = (length(egos) == n), K = K)
}

# Pooled profile from a (subset of) sources: Q_k, mean shell shares, eligible-ego shares, P_k
profile_stats <- function(prof, rows = NULL) {
  SZ <- prof$SZ; DS <- prof$DS; de <- prof$deg[prof$egos]
  if (!is.null(rows)) { SZ <- SZ[rows, , drop = FALSE]; DS <- DS[rows, , drop = FALSE]; de <- de[rows] }
  mu1 <- mean(prof$deg); mu2 <- mean(prof$deg^2); K <- ncol(SZ)
  B <- colSums(DS) / pmax(colSums(SZ), 1e-300); B[1] <- mu2 / mu1          # B_1 exact (Theorem 1(b))
  B[colSums(SZ) == 0] <- NA
  elig <- colMeans(SZ > 0)
  P <- sapply(seq_len(K), function(k) { ok <- SZ[, k] > 0; if (!any(ok)) NA else mean(DS[ok, k] / SZ[ok, k] > de[ok]) })
  list(Q = B / mu1, shell_share = colMeans(SZ) / prof$n, elig = elig, P = P, mu1 = mu1, mu2 = mu2, nu = mu2 / mu1 - 1)
}

# Reversal distances from an enrichment profile Q (Q[1] = Q_1)
k_rev_discrete <- function(Q) { j <- which(Q[-1] < 1)[1]; if (is.na(j)) NA_real_ else j + 1 }
k_rev_cont <- function(Q, method = c("linear", "log")) {
  method <- match.arg(method); k <- k_rev_discrete(Q); if (is.na(k)) return(NA_real_)
  if (method == "log") { a <- log(Q[k - 1]); b <- log(Q[k]) } else { a <- Q[k - 1] - 1; b <- Q[k] - 1 }   # crossing of Q = 1  <=>  a, b measured from 0
  (k - 1) + a / (a - b)       # linear in k (or in log Q) between k-1 (Q > 1) and k (Q < 1)
}
# Distance at which P_k crosses 0.5 (interpolated), first k with P_k < 0.5
k_P50 <- function(P) { j <- which(P < 0.5)[1]; if (is.na(j)) NA_real_ else if (j == 1) 1 else (j - 1) + (P[j - 1] - 0.5) / (P[j - 1] - P[j]) }

# Bootstrap SE of the continuous crossing (resample sources)
boot_ktil_se <- function(prof, B = 200, method = "linear") {
  ns <- nrow(prof$SZ); if (prof$exact && ns > 4000) return(NA_real_)
  v <- replicate(B, { r <- sample.int(ns, replace = TRUE); k_rev_cont(profile_stats(prof, r)$Q, method) }); sd(v, na.rm = TRUE) }

# Degree-class mean-field model (B, Sec 2.4). Input: integer degree sequence. Output: Q_k (k=1..K) and, per source degree s, the shell mean m_k(s).
mf_model <- function(deg, K = 30) {
  n <- length(deg); tab <- table(deg); dv <- as.numeric(names(tab)); nd <- as.numeric(tab); S <- sum(deg); mu1 <- S / n; U <- length(dv)
  size <- dsum <- matrix(0, U, K)
  for (s in seq_len(U)) { Fk <- dv[s]; cum <- 0; rprev <- rep(0, U)
    for (k in seq_len(K)) { cum <- cum + Fk; r <- 1 - exp(-dv * cum / S); newp <- pmax(r - rprev, 0)
      size[s, k] <- sum(nd * newp); dsum[s, k] <- sum(nd * newp * dv); Fk <- sum(nd * newp * (dv - 1)); rprev <- r; if (Fk < 1e-12) break } }
  Q <- colSums(nd * dsum) / pmax(colSums(nd * size), 1e-300) / mu1
  m <- dsum / pmax(size, 1e-300); m[size < 1e-9] <- NA                      # expected shell mean degree for a degree-s ego
  P <- sapply(seq_len(K), function(k) { ok <- !is.na(m[, k]); if (!any(ok)) NA else sum(nd[ok & m[, k] > dv]) / sum(nd[ok]) })
  thr <- sapply(seq_len(K), function(k) { ok <- !is.na(m[, k]) & m[, k] > dv; if (any(ok)) max(dv[ok]) else NA })   # largest ego degree that still has the paradox
  list(Q = Q, P = P, threshold = thr, dv = dv, nd = nd, m = m)
}
mf_krev <- function(deg, K = 30, method = "linear") { m <- mf_model(deg, K); k_rev_cont(m$Q[seq_len(max(which(!is.na(m$Q))))], method) }

# Per-graph record used by Stage 1 / Stage 2 (one row, fixed set of columns padded to KMAX)
graph_record <- function(g, nsrc = 2000, exact_n = 4000, kmax_out = 16, boot = 200, mf = TRUE, seed = 1) {
  set.seed(seed); n <- vcount(g); deg <- degree(g)
  egos <- if (n <= exact_n) NULL else sample.int(n, nsrc)
  prof <- shell_profile(g, egos, chunk = if (n > 1e5) 40 else 200); st <- profile_stats(prof)
  pad <- function(x, nm) { v <- rep(NA_real_, kmax_out); v[seq_len(min(length(x), kmax_out))] <- x[seq_len(min(length(x), kmax_out))]; setNames(v, paste0(nm, seq_len(kmax_out))) }
  kd <- k_rev_discrete(st$Q); kl <- k_rev_cont(st$Q, "linear"); kg <- k_rev_cont(st$Q, "log")
  se <- if (is.na(kl) || prof$exact || boot <= 0) NA_real_ else boot_ktil_se(prof, boot)
  mfo <- if (mf) mf_model(deg, K = 30) else NULL
  kmf <- if (mf) { Qm <- mfo$Q[seq_len(max(which(!is.na(mfo$Q) & mfo$Q > 0)))]; k_rev_cont(Qm, "linear") } else NA_real_
  cbind(data.frame(n = n, mean_deg = st$mu1, cv = sd(deg) / st$mu1, nu = st$nu, z = log(n) / log(st$nu),
                   clustering = transitivity(g, type = "global"), assort = suppressWarnings(assortativity_degree(g, directed = FALSE)),
                   K = prof$K, exact = prof$exact, krev = kd, ktil = kl, ktil_log = kg, ktil_se = se, ktil_mf = kmf, kP50 = k_P50(st$P)),
        t(pad(st$Q, "Q")), t(pad(st$P, "P")), t(pad(st$shell_share, "c")), t(pad(st$elig, "e")))
}
