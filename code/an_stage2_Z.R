# Stage 2 analysis (H5, H6, side effects). Usage: Rscript an_stage2_Z.R [file]
suppressMessages(library(dplyr)); options(width = 200, digits = 3)
f <- commandArgs(trailingOnly = TRUE); f <- if (length(f)) f[1] else "../results/stage2_full.csv"; r <- read.csv(f)
r$t_cell <- ifelse(is.na(r$t_target), "start", sprintf("%.2f", r$t_target)); cat("graphs:", nrow(r), " bases:", paste(unique(r$base), collapse = ","), "\n")
cat("\n== achieved vs target (max |r_ach - r_target|; clustering achieved by cell) ==\n"); cat("max |r_ach - r_target|:", round(max(abs(r$r_ach - r$r_target)), 4), "\n")
print(as.data.frame(r %>% group_by(base, t_cell) %>% summarise(t_ach = mean(t_ach), t_start = mean(t_start), max_dev_t = max(abs(t_ach - ifelse(is.na(t_target), t_start, t_target))), .groups = "drop")))
cat("\n== cell means: continuous crossing ktil by assortativity target (rows) and clustering cell (columns), per base ==\n")
for (b in unique(r$base)) { cat("\nbase", b, "\n"); w <- r[r$base == b, ] %>% group_by(r_target, t_cell) %>% summarise(ktil = mean(ktil), .groups = "drop") %>% tidyr::pivot_wider(names_from = t_cell, values_from = ktil); print(as.data.frame(w)) }
cat("\n== side effects by clustering cell (mean over assortativity targets): modularity, communities, mean distance ==\n")
print(as.data.frame(r %>% group_by(base, t_cell) %>% summarise(t_ach = mean(t_ach), modularity = mean(modularity), n_comm = mean(n_comm), mean_dist = mean(mean_dist), .groups = "drop")))
cat("\n== regressions of ktil (per base) ==\n")
for (b in unique(r$base)) { d <- r[r$base == b, ]
  m1 <- lm(ktil ~ r_ach + t_ach, d); m2 <- lm(ktil ~ r_ach + t_ach + modularity, d); m3 <- lm(ktil ~ r_ach * t_ach, d)
  cat("\nbase", b, " (N =", nrow(d), ")\n  ktil ~ r + T:     r %+.3f (se %.3f)  T %+.3f (se %.3f)  R2 %.3f\n" |> sprintf(coef(m1)["r_ach"], coef(summary(m1))["r_ach", 2], coef(m1)["t_ach"], coef(summary(m1))["t_ach", 2], summary(m1)$r.squared))
  cat(sprintf("  + modularity:      r %+.3f  T %+.3f  modularity %+.3f (se %.3f)  R2 %.3f\n", coef(m2)["r_ach"], coef(m2)["t_ach"], coef(m2)["modularity"], coef(summary(m2))["modularity", 2], summary(m2)$r.squared))
  cat(sprintf("  r x T interaction: %+.3f (se %.3f)\n", coef(m3)["r_ach:t_ach"], coef(summary(m3))["r_ach:t_ach", 2])) }
cat("\n== departure from the mean-field crossing (ktil - ktil_mf) by cell ==\n"); print(as.data.frame(r %>% group_by(base, r_target, t_cell) %>% summarise(dev = mean(ktil - ktil_mf), .groups = "drop") %>% tidyr::pivot_wider(names_from = t_cell, values_from = dev)))
