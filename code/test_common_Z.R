source("common_Z.R")
chk <- function(name, ok) cat(sprintf("%-70s %s\n", name, if (isTRUE(ok)) "PASS" else "FAIL"))
# 1. synthetic crossings
chk("k_rev_cont linear: Q=(1.5,1.2,0.8) -> 2.5", abs(k_rev_cont(c(1.5,1.2,0.8)) - 2.5) < 1e-12)
chk("k_rev_discrete: Q=(1.5,1.2,0.8) -> 3", k_rev_discrete(c(1.5,1.2,0.8)) == 3)
chk("k_rev: no crossing -> NA", is.na(k_rev_cont(c(1.5,1.2,1.1))))
chk("k_rev_cont log: Q=(2,1.5,0.5): 2 + log(1.5)/(log1.5-log0.5)", abs(k_rev_cont(c(2,1.5,0.5),"log") - (2 + log(1.5)/(log(1.5)-log(0.5)))) < 1e-12)
chk("k_P50 interpolation: P=(0.8,0.6,0.4) -> 2.5", abs(k_P50(c(0.8,0.6,0.4)) - 2.5) < 1e-12)
# 2. regular graph: Q_k = 1 for all k (no reversal)
g <- make_ring(30); st <- profile_stats(shell_profile(g)); chk("ring: all Q_k = 1", all(abs(st$Q[!is.na(st$Q)] - 1) < 1e-12)); chk("ring: k_rev is NA", is.na(k_rev_discrete(st$Q)))
# 3. Theorem 1 identities on several graphs (exact profile)
check_thm <- function(g, label) {
  g <- lcc(g); n <- vcount(g); p <- shell_profile(g); st <- profile_stats(p); d <- degree(g)
  B <- st$Q * st$mu1; En <- colMeans(p$SZ)
  chk(paste(label, ": B_1 = mu2/mu1"), abs(B[1] - mean(d^2)/mean(d)) < 1e-10)
  chk(paste(label, ": direct B_1 from shell contents equals mu2/mu1"), abs(sum(p$DS[,1])/sum(p$SZ[,1]) - mean(d^2)/mean(d)) < 1e-10)
  chk(paste(label, ": conservation sum_k E[n_k](B_k - mu1) = 0"), abs(sum(En * (B - st$mu1))) < 1e-8)
  # identity (a): B_k = E[d n_k]/E[n_k], with shells recomputed from the distance matrix
  D <- distances(g); ok <- TRUE; for (k in seq_len(p$K)) { nk <- rowSums(D == k); ok <- ok && abs(B[k] - sum(d*nk)/sum(nk)) < 1e-10 }
  chk(paste(label, ": identity (a) B_k = E[d n_k]/E[n_k] for all k"), ok)
  chk(paste(label, ": some Q_k < 1 for k>=2 (non-regular)"), !is.na(k_rev_discrete(st$Q)))
}
set.seed(1); check_thm(make_star(12), "star"); check_thm(add_edges(make_ring(8), c(1,5)), "ring+chord")
check_thm(sample_gnp(200, 0.03), "ER(200,.03)"); check_thm(sample_pa(150, m = 2, directed = FALSE), "BA(150,2)"); check_thm(make_tree(40, 3), "tree(40,3)")
# 4. P_k on a star: leaves' shell-1 mean degree > own; hub's shell-1 mean < own
g <- make_star(10); p <- shell_profile(g); st <- profile_stats(p); chk("star: P_1 = 9/10 (all leaves yes, hub no)", abs(st$P[1] - 0.9) < 1e-12)
# 5. sampled-vs-exact agreement
set.seed(2); g <- lcc(sample_gnp(1500, 0.02)); ex <- profile_stats(shell_profile(g)); sm <- profile_stats(shell_profile(g, sample.int(1500, 400)))
chk("sampled (400 egos) Q_2 within 0.05 of exact", abs(ex$Q[2] - sm$Q[2]) < 0.05); chk("sampled Q_1 equals exact Q_1 (B_1 set exactly)", abs(ex$Q[1] - sm$Q[1]) < 1e-12)
# 6. mean-field: Q_1 = mu2/mu1^2, threshold monotone, runs on a sequence
set.seed(3); d <- degree(lcc(sample_gnp(800, 0.03))); mf <- mf_model(d, K = 12); chk("MF: Q_1 = mu2/mu1^2", abs(mf$Q[1] - mean(d^2)/mean(d)^2) < 0.02)
thr <- mf$threshold[!is.na(mf$threshold)]; chk("MF: paradox threshold degree non-increasing in k", all(diff(thr) <= 1e-9))
# 7. graph_record end to end
set.seed(4); g <- lcc(sample_gnp(600, 0.03)); r <- graph_record(g, boot = 20); chk("graph_record: row with ktil, ktil_mf finite", is.finite(r$ktil) && is.finite(r$ktil_mf))
cat("record columns:", ncol(r), "\n"); print(round(unlist(r[1, c("n","mean_deg","cv","nu","z","krev","ktil","ktil_mf","kP50")]), 3))
