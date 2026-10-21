# Does the clustering-dependent sign of assortativity run through detected modularity?
suppressMessages(library(dplyr))
d <- read.csv("../results/stage2_full.csv"); d$ass <- d$assort
for (b in c("A","B")) { s <- d[d$base == b,]
  m0 <- lm(ktil ~ ass * clustering, s); m1 <- lm(ktil ~ ass * clustering + modularity, s); m2 <- lm(ktil ~ ass + clustering + modularity, s)
  cat(sprintf("\nBase %s (n=%d)\n", b, nrow(s)))
  cat(sprintf("  ass x clustering only      : ass %.2f, interaction %.2f, R2 %.3f\n", coef(m0)["ass"], coef(m0)["ass:clustering"], summary(m0)$r.squared))
  cat(sprintf("  + detected modularity      : ass %.2f, interaction %.2f, modularity %.2f, R2 %.3f\n", coef(m1)["ass"], coef(m1)["ass:clustering"], coef(m1)["modularity"], summary(m1)$r.squared))
  cat(sprintf("  no interaction, + modularity: ass %.2f (se %.2f), clustering %.2f, modularity %.2f, R2 %.3f\n", coef(m2)["ass"], summary(m2)$coef["ass",2], coef(m2)["clustering"], coef(m2)["modularity"], summary(m2)$r.squared))
  # effect of assortativity (+0.3 vs -0.3) at each clustering target, with and without adjusting for modularity
  for (t in c(0.1, 0.2, 0.3)) { u <- s[!is.na(s$t_target) & abs(s$t_target - t) < 1e-6 & abs(s$r_target) > 0.29,]; u$hi <- u$r_target > 0
    raw <- diff(tapply(u$ktil, u$hi, mean)); adj <- coef(lm(ktil ~ hi + modularity, u))["hiTRUE"]; dm <- diff(tapply(u$modularity, u$hi, mean))
    cat(sprintf("  clustering %.1f: assortativity +0.3 vs -0.3: raw %+.2f ; adjusted for modularity %+.2f ; modularity difference %+.2f\n", t, raw, adj, dm)) } }
# the same in Facebook
fb <- read.csv("../results/fb100_Z.csv"); fb$off <- fb$ktil - fb$ktil_mf
mf0 <- lm(off ~ assort + clustering, fb); mf1 <- lm(off ~ assort + clustering + modularity, fb); mf2 <- lm(off ~ assort * clustering + modularity, fb)
cat(sprintf("\nFacebook offset: assort %.2f -> %.2f after adding modularity; interaction with clustering (with modularity) %.2f (p=%.2f)\n", coef(mf0)["assort"], coef(mf1)["assort"], coef(mf2)["assort:clustering"], summary(mf2)$coef["assort:clustering",4]))
cat("Facebook range of assort:", round(range(fb$assort),2), " sd", round(sd(fb$assort),3), "; experiments span -0.3..0.3\n")
