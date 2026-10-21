# WP6: effect of assortativity at fixed planted partition (block-pair counts fixed), fixed degree sequence, fixed clustering.
suppressMessages(library(dplyr))
d <- read.csv("../results/stage6_full.csv"); d$cl <- ifelse(is.na(d$t_target), "natural", "0.20")
cat("graphs", nrow(d), "| max |r_ach - r_target|", signif(max(abs(d$r_ach - d$r_target)),2), "| max |planted modularity change|", signif(max(abs(d$planted_mod - d$planted_mod_start)),2), "| max |clustering change| in t-targeted cells", signif(max(abs(d$t_ach - d$t_start)[!is.na(d$t_target)] - 0),2), "\n")
cat("clustering achieved (natural cells):", round(range(d$t_ach[is.na(d$t_target)]),3), "; targeted:", round(range(d$t_ach[!is.na(d$t_target)]),3), "\n")
ci <- function(x, y) { dm <- mean(y) - mean(x); se <- sqrt(var(x)/length(x) + var(y)/length(y)); sprintf("%+.3f [%+.3f, %+.3f]", dm, dm - 1.96*se, dm + 1.96*se) }
cat("\nEffect of assortativity (+0.3 vs -0.3) on the continuous reversal distance:\n")
res <- d %>% group_by(base, mix, cl) %>% summarise(effect = ci(ktil[r_target == -0.3], ktil[r_target == 0.3]), d_louvain = round(mean(louvain_mod[r_target == 0.3]) - mean(louvain_mod[r_target == -0.3]), 3), d_meandist = round(mean(mean_dist[r_target == 0.3]) - mean(mean_dist[r_target == -0.3]), 3), .groups = "drop")
print(as.data.frame(res), right = FALSE)
cat("\nPooled over community strengths, by clustering:\n")
print(as.data.frame(d %>% group_by(base, cl) %>% summarise(effect = ci(ktil[r_target == -0.3], ktil[r_target == 0.3]), .groups = "drop")), right = FALSE)
cat("\nInteraction with clustering (r x natural/0.20), within base, fixed effects for mix:\n")
for (b in c("A","B")) { s <- d[d$base == b & d$r_target != 0,]; s$hi <- as.numeric(s$r_target > 0); s$t20 <- as.numeric(!is.na(s$t_target)); m <- lm(ktil ~ factor(mix) + hi * t20, s); print(round(summary(m)$coef[c("hi","hi:t20"),], 3)) }
write.csv(res, "../results/stage6_effects.csv", row.names = FALSE)
