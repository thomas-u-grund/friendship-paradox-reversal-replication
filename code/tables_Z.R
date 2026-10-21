suppressMessages(library(dplyr)); options(digits = 3)
s1 <- read.csv("../results/stage1_full.csv"); fb <- read.csv("../results/fb100_Z.csv"); star <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
cat("%% TABLE 1: Stage 1, NB, CV 1.0: mean crossing (simulation / mean-field) by n and mean degree\n")
d <- s1[s1$family == "NB" & s1$cv_target == 1.0, ] %>% group_by(n_target, mu_target) %>% summarise(sim = mean(ktil), mf = mean(ktil_mf), .groups = "drop")
for (n0 in sort(unique(d$n_target))) { row <- sapply(c(6, 12, 25, 50, 100), function(m) { x <- d[d$n_target == n0 & d$mu_target == m, ]; if (nrow(x)) sprintf("%.2f / %.2f", x$sim, x$mf) else "---" }); cat(sprintf("%s & %s\\\\\n", format(n0, big.mark = "{,}"), paste(row, collapse = " & "))) }
cat("\n%% TABLE 4: Facebook regressions\n")
fb$lnN <- log(fb$n); fb$mu10 <- fb$mean_deg / 10; fb$gap <- fb$ktil - fb$ktil_mf
fm <- list(f1 = lm(ktil ~ cv + assort + clustering + mu10 + lnN, fb), f2 = lm(ktil ~ cv + assort + clustering + mu10 + lnN + modularity, fb), g = lm(gap ~ cv + assort + clustering + mu10 + lnN + modularity, fb))
nm <- c("(Intercept)" = "Intercept", cv = "Degree CV", assort = "Assortativity", clustering = "Clustering", mu10 = "Mean degree (per 10)", lnN = "log(N)", modularity = "Modularity")
for (v in names(nm)) { cat(nm[v]); for (m in fm) { cf <- coef(summary(m)); if (v %in% rownames(cf)) cat(sprintf(" & %.3f%s (%.3f)", cf[v, 1], star(cf[v, 4]), cf[v, 2])) else cat(" & ---") }; cat(" \\\\\n") }
cat("$R^2$"); for (m in fm) cat(sprintf(" & %.3f", summary(m)$r.squared)); cat(" \\\\\n")
cat(sprintf("\ngap model with only assortativity and modularity: R2 = %.3f\n", summary(lm(gap ~ assort + modularity, fb))$r.squared))
cat(sprintf("VIFs (f2): %s\n", paste(sapply(c("cv", "assort", "clustering", "mu10", "lnN", "modularity"), function(v) { r2 <- summary(lm(as.formula(paste(v, "~ cv + assort + clustering + mu10 + lnN + modularity")), fb[, c("cv", "assort", "clustering", "mu10", "lnN", "modularity")]))$r.squared; sprintf("%s %.1f", v, 1/(1-r2)) }), collapse = "; ")))
