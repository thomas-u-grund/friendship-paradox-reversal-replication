# Analysis of Stage 1 per ANALYSIS_PLAN.md (H1-H4, H8; predictions P1-P4). Usage: Rscript an_stage1_Z.R [file]
suppressMessages(library(dplyr)); options(width = 200, digits = 3)
f <- commandArgs(trailingOnly = TRUE); f <- if (length(f)) f[1] else "../results/stage1_full.csv"; r <- read.csv(f)
r$logmu <- log(r$mean_deg); r$logn <- log(r$n); cat("graphs:", nrow(r), " with crossing:", sum(is.finite(r$ktil)), " exact:", sum(r$exact), "\n")
cat("\n== discrete reversal distance (all graphs) ==\n"); print(table(r$krev, useNA = "ifany"))
cat("\n== P4: mean-field vs simulation (continuous crossing) ==\n"); ok <- is.finite(r$ktil) & is.finite(r$ktil_mf)
cat(sprintf("N=%d  corr %.3f  RMSE %.3f  mean(MF - sim) %.3f  max |diff| %.3f\n", sum(ok), cor(r$ktil[ok], r$ktil_mf[ok]), sqrt(mean((r$ktil_mf[ok]-r$ktil[ok])^2)), mean(r$ktil_mf[ok]-r$ktil[ok]), max(abs(r$ktil_mf[ok]-r$ktil[ok]))))
cat("\n== P1/H3: within-family regression ktil ~ z (z = log n / log nu); slope with 95% CI ==\n")
for (fam in unique(r$family)) { d <- r[r$family == fam & is.finite(r$ktil), ]; if (nrow(d) < 5) next; m <- lm(ktil ~ z, d); ci <- confint(m)["z", ]
  cat(sprintf("%-3s N=%4d slope %.3f [%.3f, %.3f]  intercept %.3f  resid SD %.3f  R2 %.3f\n", fam, nrow(d), coef(m)["z"], ci[1], ci[2], coef(m)["(Intercept)"], sigma(m), summary(m)$r.squared)) }
cat("\n  by NB cv_target (realised cv in brackets):\n")
for (cvt in sort(unique(r$cv_target[r$family == "NB"]))) { d <- r[r$family == "NB" & r$cv_target == cvt & is.finite(r$ktil), ]; if (nrow(d) < 5) next; m <- lm(ktil ~ z, d); ci <- confint(m)["z", ]
  cat(sprintf("  NB cv %.1f [%.2f] N=%4d slope %.3f [%.3f, %.3f] intercept %.3f resid SD %.3f\n", cvt, mean(d$cv), nrow(d), coef(m)["z"], ci[1], ci[2], coef(m)["(Intercept)"], sigma(m))) }
cat("\n== H1/H2: Negative-Binomial, ktil ~ log(mean degree) + log n + cv (fixed effects of cv target) ==\n")
d <- r[r$family == "NB" & is.finite(r$ktil), ]; m <- lm(ktil ~ logmu + logn + factor(cv_target), d); print(round(coef(summary(m))[c("logmu", "logn"), ], 3)); cat("R2", round(summary(m)$r.squared, 3), "\n")
cat("\n== H4: cell means of ktil by n, mean degree, cv target; change from lowest to highest cv ==\n")
cm <- d %>% group_by(n_target, mu_target, cv_target) %>% summarise(k = n(), ktil = mean(ktil), se = sd(ktil) / sqrt(n()), mf = mean(ktil_mf), .groups = "drop")
w <- cm %>% select(n_target, mu_target, cv_target, ktil) %>% tidyr::pivot_wider(names_from = cv_target, values_from = ktil, names_prefix = "cv"); cvc <- grep("^cv", names(w)); w$change <- w[[cvc[length(cvc)]]] - w[[cvc[1]]]; print(as.data.frame(w))
cat("\n  sign of change by mean degree (share of n-levels with positive change):\n"); print(w %>% group_by(mu_target) %>% summarise(levels = n(), share_positive = mean(change > 0, na.rm = TRUE), mean_change = mean(change, na.rm = TRUE)))
cat("\n== P3/H8: kP50 - ktil ==\n"); ok <- is.finite(r$kP50) & is.finite(r$ktil); cat(sprintf("N=%d  mean %.3f  SD %.3f  corr(kP50, ktil) %.3f\n", sum(ok), mean(r$kP50[ok]-r$ktil[ok]), sd(r$kP50[ok]-r$ktil[ok]), cor(r$kP50[ok], r$ktil[ok])))
cat("\n== interpolation robustness: linear vs log ==\n"); ok <- is.finite(r$ktil) & is.finite(r$ktil_log); cat(sprintf("mean diff %.3f, max |diff| %.3f\n", mean(r$ktil[ok]-r$ktil_log[ok]), max(abs(r$ktil[ok]-r$ktil_log[ok]))))
cat("\n== sampling noise (sampled graphs only): bootstrap SE of ktil ==\n"); s <- r[!r$exact & is.finite(r$ktil_se), ]; if (nrow(s)) cat(sprintf("N=%d  mean SE %.3f  max SE %.3f\n", nrow(s), mean(s$ktil_se), max(s$ktil_se)))
