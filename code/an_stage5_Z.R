suppressMessages(library(dplyr)); options(width = 200, digits = 3); r <- read.csv("../results/stage5_full.csv"); r$t_cell <- ifelse(is.na(r$t_target), "start", sprintf("%.2f", r$t_target))
cmp <- function(a, b) { t <- t.test(a, b); c(diff = mean(a) - mean(b), lo = t$conf.int[1], hi = t$conf.int[2]) }
cat("graphs:", nrow(r), "\n\n== design check: planted modularity unchanged by rewiring (max |change|), assortativity unchanged (max |change|), clustering target error ==\n")
cat(sprintf("%.2e  %.2e  %.3f\n", max(abs(r$planted_mod - r$planted_mod_start)), max(abs(r$r_ach - r$r_start)), max(abs(r$t_ach - ifelse(is.na(r$t_target), r$t_start, r$t_target)))))
cat("\n== cell means ==\n"); print(as.data.frame(r %>% group_by(base, mix, t_cell) %>% summarise(k = n(), planted_Q = mean(planted_mod), louvain_Q = mean(louvain_mod), T = mean(t_ach), dist = mean(mean_dist), ktil = mean(ktil), .groups = "drop")))
cat("\n== clustering effect at FIXED planted community structure: highest clustering cell (0.30) minus lowest cell, by base and mixing; with Louvain modularity change ==\n")
for (b in c("A", "B")) for (m in sort(unique(r$mix))) { d <- r[r$base == b & r$mix == m, ]; lo <- d[d$t_ach <= quantile(d$t_ach, 0.25), ]; hi <- d[d$t_ach >= quantile(d$t_ach, 0.75), ]; x <- cmp(hi$ktil, lo$ktil)
  cat(sprintf("base %s mix %.2f (planted Q %.2f): T %.2f -> %.2f; ktil %+.3f [%+.3f, %+.3f]; Louvain Q %.2f -> %.2f\n", b, m, mean(d$planted_mod), mean(lo$t_ach), mean(hi$t_ach), x[1], x[2], x[3], mean(lo$louvain_mod), mean(hi$louvain_mod))) }
cat("\n== community effect at FIXED clustering (start cell: no triangle targeting): ktil by mixing ==\n"); print(as.data.frame(r %>% filter(t_cell == "start") %>% group_by(base, mix) %>% summarise(planted_Q = mean(planted_mod), T = mean(t_ach), ktil = mean(ktil), .groups = "drop")))
cat("\n== regressions per base ==\n")
for (b in c("A", "B")) { d <- r[r$base == b, ]
  m1 <- lm(ktil ~ t_ach + planted_mod, d); m2 <- lm(ktil ~ t_ach + planted_mod + louvain_mod, d); m3 <- lm(ktil ~ t_ach * planted_mod, d); m4 <- lm(ktil ~ t_ach + planted_mod + mean_dist, d)
  cat(sprintf("\nbase %s (N=%d)\n  ktil ~ T + planted Q:           T %+.3f (se %.3f)  Q %+.3f (se %.3f)  R2 %.3f\n", b, nrow(d), coef(m1)["t_ach"], coef(summary(m1))["t_ach", 2], coef(m1)["planted_mod"], coef(summary(m1))["planted_mod", 2], summary(m1)$r.squared))
  cat(sprintf("  + Louvain Q:                   T %+.3f  planted Q %+.3f  Louvain Q %+.3f (se %.3f)  R2 %.3f\n", coef(m2)["t_ach"], coef(m2)["planted_mod"], coef(m2)["louvain_mod"], coef(summary(m2))["louvain_mod", 2], summary(m2)$r.squared))
  cat(sprintf("  T x planted Q interaction:     %+.3f (se %.3f)\n  + mean distance:               T %+.3f  planted Q %+.3f  distance %+.3f  R2 %.3f\n", coef(m3)["t_ach:planted_mod"], coef(summary(m3))["t_ach:planted_mod", 2], coef(m4)["t_ach"], coef(m4)["planted_mod"], coef(m4)["mean_dist"], summary(m4)$r.squared))
  cat(sprintf("  cor(T, planted Q) %.3f; cor(T, Louvain Q) %.3f; cor(T, mean dist) %.3f\n", cor(d$t_ach, d$planted_mod), cor(d$t_ach, d$louvain_mod), cor(d$t_ach, d$mean_dist))) }
