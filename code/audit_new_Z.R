# Audit 2: numbers added after the first audit (coverage diagnostic, Stage 6, Facebook clustering regime, leave-one-out, slopes, proof audits).
suppressMessages(library(dplyr)); tex <- paste(readLines("../manuscript.tex"), collapse = "\n"); tex <- gsub("\\s+", " ", gsub("\\$", "", tex)); tex <- gsub("\\{,\\}", ",", tex)
fails <- 0; n <- 0
chk <- function(label, ...) { for (s in c(...)) { n <<- n + 1; if (!grepl(s, tex, fixed = TRUE)) { fails <<- fails + 1; cat("MISSING [", label, "]:", s, "\n") } } }
f2 <- function(x) sprintf("%.2f", x); f3 <- function(x) sprintf("%.3f", x)
cov_at <- function(d) { K <- 16; cs <- as.matrix(d[, paste0("c", 1:K)]); cs[is.na(cs)] <- 0; R <- t(apply(cs, 1, cumsum)); R <- cbind(0, R)
  mapply(function(i) { k <- d$ktil[i]; if (is.na(k)) return(NA_real_); lo <- floor(k); hi <- lo + 1; if (hi > K) return(NA_real_); R[i, lo + 1] + (k - lo) * (R[i, hi + 1] - R[i, lo + 1]) }, seq_len(nrow(d))) }
s1 <- read.csv("../results/stage1_full.csv"); s2 <- read.csv("../results/stage2_full.csv"); s5 <- read.csv("../results/stage5_full.csv"); s6 <- read.csv("../results/stage6_full.csv")
fb <- read.csv("../results/fb100_Z.csv"); ex <- read.csv("../results/external_Z.csv")
# ---- coverage at the crossing ----
s1$Rx <- cov_at(s1); s2$Rx <- cov_at(s2); s5$Rx <- cov_at(s5)
chk("coverage stage1", f2(mean(s1$Rx, na.rm = TRUE)), f2(sd(s1$Rx, na.rm = TRUE)), f2(sd(s1$ktil)))
for (b in c("A", "B")) { a <- s5[s5$base == b,]; z <- s2[s2$base == b,]; chk(paste("coverage designs", b), paste0(f2(mean(a$Rx, na.rm = TRUE)), " (SD ", f2(sd(a$Rx, na.rm = TRUE)), ")"), paste0(f2(mean(z$Rx, na.rm = TRUE)), " (SD ", f2(sd(z$Rx, na.rm = TRUE)), ")")) }
eff <- function(d, v) 0.1 * coef(lm(d[[v]] ~ clustering, d))[2]
ek <- c(); er <- c(); gr <- c()
for (dd in list(s5, s2)) for (b in c("A", "B")) { s <- dd[dd$base == b & !is.na(dd$Rx),]; ek <- c(ek, eff(s, "ktil")); er <- c(er, eff(s, "Rx")); s$g12 <- s$c2 / s$c1; gr <- c(gr, -eff(s, "g12") / mean(s$g12)) }
chk("coverage effects", f2(min(ek)), f2(max(ek)), f2(min(er)), f2(max(er)), paste0(round(100 * min(gr)), " to ", round(100 * max(gr)), "\\%") )
# ---- deciles ----
for (b in c("A", "B")) { s <- s5[s5$base == b,]; lo <- s[s$clustering <= quantile(s$clustering, .1),]; hi <- s[s$clustering >= quantile(s$clustering, .9),]; cat(sprintf("  info base %s deciles: low krev %s | high krev %s | mean clustering %.2f / %.2f\n", b, paste(names(table(lo$krev)), table(lo$krev), collapse = " "), paste(names(table(hi$krev)), table(hi$krev), collapse = " "), mean(lo$clustering), mean(hi$clustering))) }
chk("deciles", "mean 0.05", "mean 0.30", "15 and 33")
# ---- Stage 6 ----
d <- s6[s6$r_target != 0,]; d$hi <- as.numeric(d$r_target > 0); d$t20 <- as.numeric(!is.na(d$t_target))
chk("stage6 design", "1,080 graphs", paste0("within ", f2(ceiling(max(abs(s6$r_ach - s6$r_target)) * 100) / 100)), "Planted modularity and the degree sequence are exactly fixed")
for (b in c("A", "B")) { s <- d[d$base == b,]; m <- lm(ktil ~ factor(mix) + hi * t20, s); m1 <- lm(ktil ~ factor(mix) + hi * t20 + louvain_mod, s); m2 <- lm(ktil ~ factor(mix) + hi * t20 + mean_dist, s)
  cat(sprintf("  info stage6 base %s: interaction %.3f (adj. modularity %.3f; adj. mean distance %.3f)\n", b, coef(m)["hi:t20"], coef(m1)["hi:t20"], coef(m2)["hi:t20"]))
  chk(paste("stage6 interaction", b), f2(abs(coef(m)["hi:t20"])), f2(abs(coef(m1)["hi:t20"])), f2(abs(coef(m2)["hi:t20"]))) }
eff6 <- read.csv("../results/stage6_effects.csv"); cat("  info stage6 cell effects at 0.20:", gsub("\\[.*", "", eff6$effect[eff6$cl == "0.20"]), "\n")
chk("stage6 contrasts", "+0.135", "-0.028", "+0.209", "+0.026") 
# ---- Facebook regime, leave-one-out, slopes ----
chk("FB clustering", paste0(f2(min(fb$clustering)), " and ", f2(max(fb$clustering))), paste0(round(100 * mean(fb$clustering >= .1 & fb$clustering <= .2)), "\\% between"))
fb$off <- fb$ktil - fb$ktil_mf; fb$cl3 <- cut(fb$clustering, quantile(fb$clustering, 0:3 / 3), include.lowest = TRUE)
sl <- sapply(levels(fb$cl3), function(l) coef(lm(off ~ assort, fb[fb$cl3 == l,]))[2]); chk("FB tercile slopes", f2(min(sl)), f2(max(sl)))
chk("FB interaction p", sub("^0", "", f2(summary(lm(off ~ assort * clustering + modularity, fb))$coef["assort:clustering", 4])))
ex$off <- ex$ktil - ex$ktil_mf; loo <- sapply(seq_len(nrow(ex)), function(i) c(cor(ex$off[-i], ex$modularity[-i], method = "spearman"), cor(ex$off[-i], ex$modularity[-i])))
chk("leave-one-out", paste0("between ", f2(min(loo[1,])), " and ", f2(ceiling(max(loo[1,]) * 100) / 100)), paste0("between ", f2(min(loo[2,])), " and ", f2(max(loo[2,]))))
sx <- function(d, v) coef(lm(d[[v]] ~ I(log(d$n) / log(d$nu))))[2]; chk("slopes", f2(sx(fb, "ktil")), f2(sx(fb, "ktil_mf")), f2(sx(ex, "ktil")), f2(sx(ex, "ktil_mf")))
# ---- proof audit and law counts quoted in the text (from the saved outputs) ----
pa <- paste(readLines("../results/proof_audit.txt"), collapse = "\n"); pb <- paste(readLines("../results/proof_audit2.txt"), collapse = "\n")
for (s in c("70/70", "6.2e-15", "1.54e-02", "9.69e-08")) { n <- n + 1; if (!grepl(s, pa, fixed = TRUE)) { fails <- fails + 1; cat("MISSING in proof_audit.txt:", s, "\n") } }
for (s in c("27/27", "Q_inf from 0.034 to 0.685", "9 with Q_inf>=1")) { n <- n + 1; if (!grepl(s, pb, fixed = TRUE)) { fails <- fails + 1; cat("MISSING in proof_audit2.txt:", s, "\n") } }
chk("proof numbers in text", "6\\times10^{-15}", "1.5\\times10^{-2}", "34 of them", "nine", "0.03 and 0.69", "400 values")
cat(sprintf("\naudit 2: %d checks, %d missing\n", n, fails))
