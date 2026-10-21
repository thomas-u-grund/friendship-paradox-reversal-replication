suppressMessages(library(dplyr)); options(width = 200, digits = 3); r <- read.csv("../results/stage2_full.csv"); r$t_cell <- ifelse(is.na(r$t_target), "start", sprintf("%.2f", r$t_target))
cmp <- function(a, b) { t <- t.test(a, b); c(diff = mean(a) - mean(b), lo = t$conf.int[1], hi = t$conf.int[2]) }
cat("== Clustering effect (T = 0.30 vs start) on ktil, by assortativity target ==\n")
for (b in c("A", "B")) for (rt in c(-0.3, 0, 0.3)) { d <- r[r$base == b & r$r_target == rt, ]; x <- cmp(d$ktil[d$t_cell == "0.30"], d$ktil[d$t_cell == "start"]); cat(sprintf("base %s r=%+.1f: %+.3f [%+.3f, %+.3f]\n", b, rt, x[1], x[2], x[3])) }
cat("\n== Assortativity effect (r=+0.3 vs r=-0.3) on ktil, by clustering cell ==\n")
for (b in c("A", "B")) for (tc in c("start", "0.10", "0.20", "0.30")) { d <- r[r$base == b & r$t_cell == tc, ]; x <- cmp(d$ktil[d$r_target == 0.3], d$ktil[d$r_target == -0.3]); cat(sprintf("base %s T=%-5s: %+.3f [%+.3f, %+.3f]\n", b, tc, x[1], x[2], x[3])) }
cat("\n== Modularity side effect (T = 0.30 vs start), and mean distance ==\n")
for (b in c("A", "B")) { d <- r[r$base == b, ]; x <- cmp(d$modularity[d$t_cell == "0.30"], d$modularity[d$t_cell == "start"]); y <- cmp(d$mean_dist[d$t_cell == "0.30"], d$mean_dist[d$t_cell == "start"]); cat(sprintf("base %s: modularity %+.3f [%+.3f, %+.3f]; mean distance %+.3f [%+.3f, %+.3f]\n", b, x[1], x[2], x[3], y[1], y[2], y[3])) }
cat("\n== Does modularity explain the clustering effect? partial R2 / mediation (per base, pooled over r) ==\n")
for (b in c("A", "B")) { d <- r[r$base == b, ]; m0 <- lm(ktil ~ t_ach + r_ach, d); m1 <- lm(ktil ~ t_ach + r_ach + modularity, d); m2 <- lm(ktil ~ modularity + r_ach, d)
  cat(sprintf("base %s: T coefficient %+.2f -> %+.2f after adding modularity; R2 %.3f -> %.3f; modularity alone (+r) R2 %.3f; cor(T, modularity) %.3f\n", b, coef(m0)["t_ach"], coef(m1)["t_ach"], summary(m0)$r.squared, summary(m1)$r.squared, summary(m2)$r.squared, cor(d$t_ach, d$modularity))) }
cat("\n== Discrete k_rev by clustering cell (share with k_rev >= 4) ==\n"); print(as.data.frame(r %>% group_by(base, t_cell) %>% summarise(share_krev_ge4 = mean(krev >= 4), mean_krev = mean(krev), .groups = "drop")))
cat("\n== Departure from mean-field vs real-network offset: cell (T=0.20, r=0..0.15) ==\n"); print(as.data.frame(r %>% filter(t_cell == "0.20", r_target %in% c(0, 0.15)) %>% group_by(base, r_target) %>% summarise(ktil_minus_mf = mean(ktil - ktil_mf), se = sd(ktil - ktil_mf) / sqrt(n()), .groups = "drop")))
