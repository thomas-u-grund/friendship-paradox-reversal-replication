# Checks prompted by external feedback (2026-10-21): clustering regime of Facebook networks; leave-one-out for the ten; shell growth rate vs clustering.
suppressMessages(library(dplyr))
fb <- read.csv("../results/fb100_Z.csv"); fb$off <- fb$ktil - fb$ktil_mf
cat("Facebook clustering: mean", round(mean(fb$clustering),3), "range", round(range(fb$clustering),3), "quantiles", round(quantile(fb$clustering,c(.1,.25,.5,.75,.9)),3), "\n")
cat("Facebook assortativity: mean", round(mean(fb$assort),3), "range", round(range(fb$assort),3), "\n")
cat("share with clustering < 0.10, 0.10-0.20, > 0.20:", round(mean(fb$clustering<.10),2), round(mean(fb$clustering>=.10&fb$clustering<=.20),2), round(mean(fb$clustering>.20),2), "\n")
fb$cl3 <- cut(fb$clustering, quantile(fb$clustering, 0:3/3), include.lowest = TRUE)
for (l in levels(fb$cl3)) { s <- fb[fb$cl3 == l,]; m <- lm(off ~ assort, s); cat(sprintf("  clustering %s (n=%d): slope of offset on assortativity %.2f (se %.2f)\n", l, nrow(s), coef(m)[2], summary(m)$coefficients[2,2])) }
m <- lm(off ~ assort * clustering, fb); print(round(summary(m)$coefficients, 3))
ex <- read.csv("../results/external_Z.csv"); ex$off <- ex$ktil - ex$ktil_mf
cat("\nTen networks: Spearman(offset, modularity) =", round(cor(ex$off, ex$modularity, method = "spearman"), 3), " Pearson =", round(cor(ex$off, ex$modularity), 3), "\n")
loo <- sapply(seq_len(nrow(ex)), function(i) c(sp = cor(ex$off[-i], ex$modularity[-i], method = "spearman"), pe = cor(ex$off[-i], ex$modularity[-i])))
cat("Leave-one-out Spearman range", round(range(loo["sp",]),3), " Pearson range", round(range(loo["pe",]),3), " (Pearson without Amazon:", round(loo["pe", which(grepl("mazon", ex$label))],3), ")\n")
print(ex[, c("label","modularity","off")])
cat("Ten: max CV", ex$label[which.max(ex$cv)], "; min mean degree", ex$label[which.min(ex$mean_deg)], "\n")
# shell growth under clustering (Stage 5 planted partitions and Stage 2 rewiring), growth ratio of shell shares c_{k+1}/c_k
for (nm in c("stage5","stage2")) { d <- read.csv(paste0("../results/", nm, "_full.csv")); d$g12 <- d$c2/d$c1; d$g23 <- d$c3/d$c2
  cat("\n", nm, ": effect of +0.1 clustering (within base) on growth ratios\n"); for (b in unique(d$base)) { s <- d[d$base==b,]
   cat(sprintf("  base %s: c2/c1 %.2f -> %.3f per 0.1 clustering (mean %.1f); c3/c2 %.3f per 0.1 (mean %.2f); c4/c3 %.3f per 0.1\n", b, 0, 0.1*coef(lm(g12~clustering,s))[2], mean(s$g12), 0.1*coef(lm(g23~clustering,s))[2], mean(s$g23, na.rm=TRUE), 0.1*coef(lm(I(c4/c3)~clustering,s))[2])) } }
