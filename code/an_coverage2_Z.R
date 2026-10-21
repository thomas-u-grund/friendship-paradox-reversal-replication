# Coverage at the crossing: R(ktil) = cumulative share of the network within distance ktil (linear interpolation), per graph.
suppressMessages(library(dplyr))
cov_at <- function(d) { K <- 16; cs <- as.matrix(d[, paste0("c", 1:K)]); cs[is.na(cs)] <- 0; R <- t(apply(cs, 1, cumsum)); R <- cbind(0, R)   # R[,j+1] = coverage within distance j (distance 0 = ego itself, ignore)
  mapply(function(i) { k <- d$ktil[i]; if (is.na(k)) return(NA_real_); lo <- floor(k); hi <- lo + 1; if (hi > K) return(NA_real_); R[i, lo + 1] + (k - lo) * (R[i, hi + 1] - R[i, lo + 1]) }, seq_len(nrow(d))) }
f <- function(x) formatC(x, format = "f", digits = 3)
for (nm in c("stage5", "stage2")) {
  d <- read.csv(paste0("../results/", nm, "_full.csv")); d$Rx <- cov_at(d); d <- d[!is.na(d$Rx) & !is.na(d$ktil), ]
  cond <- if (nm == "stage5") c("mix") else c("r_target")
  cat("\n==", nm, "(n graphs", nrow(d), ")\n")
  for (b in unique(d$base)) { s <- d[d$base == b, ]
    cat(sprintf("base %s: ktil mean %.2f sd %.2f range %.2f-%.2f | coverage-at-crossing mean %.3f sd %.3f range %.3f-%.3f | cor(ktil, R_cross) %.2f\n", b, mean(s$ktil), sd(s$ktil), min(s$ktil), max(s$ktil), mean(s$Rx), sd(s$Rx), min(s$Rx), max(s$Rx), cor(s$ktil, s$Rx)))
    m <- lm(ktil ~ clustering + assort, data = s); m2 <- lm(Rx ~ clustering + assort, data = s)
    cat(sprintf("   effect of +0.1 clustering: on ktil %.2f (R2 %.2f) ; on coverage-at-crossing %.3f (R2 %.2f)\n", 0.1 * coef(m)["clustering"], summary(m)$r.squared, 0.1 * coef(m2)["clustering"], summary(m2)$r.squared)) }
  # relative variability: how much of the spread in ktil is spread in coverage? Compare with the between-base spread
  write.csv(d[, c("base", cond, "clustering", "assort", "ktil", "Rx")], paste0("../results/coverage_cross_", nm, ".csv"), row.names = FALSE) }
# Stage 1 (degree law and size vary): coverage at crossing vs degree law
d <- read.csv("../results/stage1_full.csv"); d$Rx <- cov_at(d); d <- d[!is.na(d$Rx) & !is.na(d$ktil), ]
cat("\nStage 1: coverage at crossing mean", f(mean(d$Rx)), "sd", f(sd(d$Rx)), "| ktil sd", f(sd(d$ktil)), "\n")
cat("  by family:\n"); print(d %>% group_by(family = if ("family" %in% names(d)) family else "all") %>% summarise(n = n(), R_mean = mean(Rx), R_sd = sd(Rx), k_sd = sd(ktil)))
