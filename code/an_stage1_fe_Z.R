# Refined H3 test (post hoc; recorded in the deviation log): fix the degree law (cell mu_target x cv_target) and vary only n.
suppressMessages(library(dplyr)); options(width = 200, digits = 3); r <- read.csv("../results/stage1_full.csv")
r$cell <- paste(r$family, r$mu_target, r$cv_target)
cat("== slope of ktil on z with degree-law fixed effects (cell intercepts), by family ==\n")
for (fam in c("NB", "LN", "PL")) { d <- r[r$family == fam, ]; m <- lm(ktil ~ z + factor(cell), d); ci <- confint(m)["z", ]; cat(sprintf("%-2s N=%d cells=%d slope %.3f [%.3f, %.3f] resid SD %.3f\n", fam, nrow(d), length(unique(d$cell)), coef(m)["z"], ci[1], ci[2], sigma(m))) }
cat("\n== same, mean-field values on the same graphs ==\n")
for (fam in c("NB", "LN", "PL")) { d <- r[r$family == fam, ]; m <- lm(ktil_mf ~ z + factor(cell), d); ci <- confint(m)["z", ]; cat(sprintf("%-2s MF slope %.3f [%.3f, %.3f]\n", fam, coef(m)["z"], ci[1], ci[2])) }
cat("\n== per-cell slopes (NB): ktil on log n / log nu within (mu, cv) cell; summary over cells ==\n")
cs <- r[r$family == "NB", ] %>% group_by(mu_target, cv_target) %>% filter(n_distinct(n_target) >= 4) %>% group_modify(~{ m <- lm(ktil ~ z, .x); data.frame(levels = n_distinct(.x$n_target), slope = coef(m)["z"], se = coef(summary(m))["z", 2]) }) %>% ungroup()
print(as.data.frame(cs)); cat(sprintf("\nNB per-cell slopes: mean %.3f, SD %.3f, min %.3f, max %.3f (cells = %d)\n", mean(cs$slope), sd(cs$slope), min(cs$slope), max(cs$slope), nrow(cs)))
cat("\n== per-cell slope vs log n only (ktil ~ log n) / (1/log nu): equivalent slope on log n ==\n")
cl <- r[r$family == "NB", ] %>% group_by(mu_target, cv_target) %>% filter(n_distinct(n_target) >= 4) %>% group_modify(~{ m <- lm(ktil ~ log(n), .x); data.frame(slope_logn = coef(m)["log(n)"], nu = mean(.x$nu), predicted = 1 / log(mean(.x$nu))) }) %>% ungroup(); cl$ratio <- cl$slope_logn / cl$predicted; cat(sprintf("slope_logn / (1/log nu): mean %.3f, range %.3f-%.3f\n", mean(cl$ratio), min(cl$ratio), max(cl$ratio)))
cat("\n== intercept c = ktil - z by cell (NB), mean over graphs ==\n"); r$c <- r$ktil - r$z; print(as.data.frame(r[r$family=="NB",] %>% group_by(mu_target, cv_target) %>% summarise(c = mean(c), sd_c = sd(c), .groups = "drop") %>% tidyr::pivot_wider(names_from = cv_target, values_from = c, names_prefix = "cv")))
