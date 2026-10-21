suppressMessages(library(dplyr)); options(width = 200, digits = 3); r <- read.csv("../results/stage1_full.csv")
nbr <- r[r$family == "NB", ]
cell <- nbr %>% group_by(mu_target, cv_target, n_target) %>% summarise(ktil = mean(ktil), mf = mean(ktil_mf), nu = mean(nu), .groups = "drop")
rat <- cell %>% group_by(mu_target, cv_target) %>% filter(n() >= 4) %>% summarise(levels = n(), nu = mean(nu), sim = coef(lm(ktil ~ log(n_target)))[2] * log(mean(nu)), mf = coef(lm(mf ~ log(n_target)))[2] * log(mean(nu)), .groups = "drop")
cat("ratio = (d ktil / d log n) * log(nu)  [=1 under the scaling statement], cell means over graphs, by (mean degree, cv)\n"); print(as.data.frame(rat))
cat(sprintf("\nsimulation: mean %.3f SD %.3f range %.3f-%.3f | mean-field on same cells: mean %.3f SD %.3f range %.3f-%.3f\n", mean(rat$sim), sd(rat$sim), min(rat$sim), max(rat$sim), mean(rat$mf), sd(rat$mf), min(rat$mf), max(rat$mf)))
cat("\nby mean degree:\n"); print(rat %>% group_by(mu_target) %>% summarise(sim = mean(sim), mf = mean(mf), .groups = "drop") %>% as.data.frame())
cat("\ncorr(sim ratio, MF ratio) across cells:", round(cor(rat$sim, rat$mf), 3), "\n")
# local slope in the asymptotic direction: compare adjacent n levels
loc <- cell %>% group_by(mu_target, cv_target) %>% arrange(n_target) %>% mutate(dk = c(NA, diff(ktil)), dlogn = c(NA, diff(log(n_target))), loc = dk / dlogn * log(nu)) %>% ungroup()
cat("\nlocal ratio (adjacent n levels), by n step: sim vs MF\n"); loc$locmf <- (loc$mf - dplyr::lag(loc$mf)) / loc$dlogn * log(loc$nu)
print(loc %>% filter(!is.na(loc)) %>% group_by(n_target) %>% summarise(sim = mean(loc), mf = mean(locmf, na.rm = TRUE), .groups = "drop") %>% as.data.frame())
