# Audit: recompute the numbers quoted in main_Z.tex from the data files and check that each appears in the text.
suppressMessages(library(dplyr)); tex <- paste(readLines("../manuscript.tex"), collapse = "\n"); tex <- gsub("\\s+", " ", tex)
fails <- 0; n <- 0
chk <- function(label, ...) { for (s in c(...)) { n <<- n + 1; if (!grepl(s, tex, fixed = TRUE)) { fails <<- fails + 1; cat("MISSING [", label, "]:", s, "\n") } } }
f2 <- function(x) sprintf("%.2f", x); f3 <- function(x) sprintf("%.3f", x); big <- function(x) format(x, big.mark = "{,}")
s1 <- read.csv("../results/stage1_full.csv"); s2 <- read.csv("../results/stage2_full.csv"); s5 <- read.csv("../results/stage5_full.csv"); fb <- read.csv("../results/fb100_Z.csv"); ego <- read.csv("ego_rel_all.csv")
# ---- Stage 1 ----
tb <- table(s1$krev); chk("krev counts", paste0(big(tb["3"]), " graphs"), paste0("4 in ", big(tb["4"])), paste0("5 in ", big(tb["5"])), paste0("6 in ", big(tb["6"]))); chk("n graphs", big(nrow(s1)), paste0(big(sum(s1$exact)), " graphs"))
chk("MF vs sim", f3(cor(s1$ktil, s1$ktil_mf)), f2(sqrt(mean((s1$ktil_mf - s1$ktil)^2))), f2(mean(s1$ktil_mf - s1$ktil)), f2(max(abs(s1$ktil_mf - s1$ktil))), f3(max(s1$ktil_se, na.rm = TRUE)))
nb <- s1[s1$family == "NB", ]; nb$lm <- log(nb$mean_deg); nb$ln <- log(nb$n); m <- lm(ktil ~ lm + ln + factor(cv_target), nb); chk("regression", f2(-coef(m)["lm"]), f2(coef(m)["ln"]), f2(summary(m)$r.squared))
d <- nb[nb$cv_target == 1, ] %>% group_by(n_target, mu_target) %>% summarise(a = mean(ktil), b = mean(ktil_mf), se = sd(ktil)/sqrt(n()), .groups = "drop"); for (i in seq_len(nrow(d))) chk("Table 1", paste0(f2(d$a[i]), " / ", f2(d$b[i]))); chk("Table 1 SE", f3(max(d$se)))
cell <- nb %>% group_by(mu_target, cv_target, n_target) %>% summarise(k = mean(ktil), mf = mean(ktil_mf), nu = mean(nu), .groups = "drop")
rat <- cell %>% group_by(mu_target, cv_target) %>% filter(n() >= 4) %>% summarise(sim = coef(lm(k ~ log(n_target)))[2] * log(mean(nu)), mf = coef(lm(mf ~ log(n_target)))[2] * log(mean(nu)), .groups = "drop")
chk("ratio", f3(mean(rat$sim)), f3(sd(rat$sim)), f2(min(rat$sim)), f2(max(rat$sim)), paste0(sum(rat$sim >= 0.9 & rat$sim <= 1.1), " of 20"), f2(mean(rat$mf)), f2(min(rat$mf)), f2(max(rat$mf)), f2(cor(rat$sim, rat$mf)))
for (fam in c("NB", "LN", "PL")) { dd <- s1[s1$family == fam, ]; mm <- lm(ktil ~ z, dd); ci <- confint(mm)["z", ]; chk(paste("pooled slope", fam), f3(coef(mm)["z"]), paste0("[", f3(ci[1]), ",", f3(ci[2]), "]")) }
ch <- cell %>% select(mu_target, cv_target, n_target, k) %>% group_by(mu_target, n_target) %>% summarise(ch = k[cv_target == 1.4] - k[cv_target == 0.4], .groups = "drop")
for (mu in c(6, 12, 25, 50, 100)) { x <- ch[ch$mu_target == mu, ]; cat(sprintf("  info: change in crossing CV .4->1.4 at mean degree %d by n: %s\n", mu, paste(sprintf("%+.2f", x$ch), collapse = " "))) }
# ---- Stage 2 / WP5 ----
s2$tc <- ifelse(is.na(s2$t_target), "start", sprintf("%.2f", s2$t_target)); ci <- function(a, b) { t <- t.test(a, b); c(unname(t$estimate[1] - t$estimate[2]), t$conf.int) }; fmt <- function(v) paste0(f3(v[2]), ", ", f3(v[3]))
for (b in c("A", "B")) { dd <- s2[s2$base == b, ]; x <- ci(dd$ktil[dd$r_target == 0 & dd$tc == "0.30"], dd$ktil[dd$r_target == 0 & dd$tc == "start"]); chk(paste("T2 clustering", b), f3(x[1]), fmt(x))
  for (tc in c("start", "0.10", "0.20", "0.30")) { x <- ci(dd$ktil[dd$r_target == 0.3 & dd$tc == tc], dd$ktil[dd$r_target == -0.3 & dd$tc == tc]); chk(paste("T2 assort", b, tc), f3(abs(x[1])), paste0(f3(abs(x[2])), ", ", f3(abs(x[3])))) } }
for (b in c("A", "B")) { dd <- s2[s2$base == b & s2$r_target == 0, ]; chk("modularity/distance", "0.34", "0.37") }
s5$tc <- ifelse(is.na(s5$t_target), "start", sprintf("%.2f", s5$t_target))
for (b in c("A", "B")) for (m in c(0.30, 0.15)) { dd <- s5[s5$base == b & s5$mix == m, ]; lo <- dd[dd$t_ach <= quantile(dd$t_ach, .25), ]; hi <- dd[dd$t_ach >= quantile(dd$t_ach, .75), ]; x <- ci(hi$ktil, lo$ktil); chk(paste("T2 WP5 clustering", b, m), f3(x[1]), fmt(x)) }
for (b in c("A", "B")) { s <- s5[s5$base == b & s5$tc == "0.10" & s5$mix == 0.15, ]; w <- s5[s5$base == b & s5$tc == "0.10" & s5$mix == 0.80, ]; x <- ci(s$ktil, w$ktil); chk(paste("T2 WP5 community", b), f3(x[1]), fmt(x)) }
chk("WP5 mean distance corr", f3(cor(s5$ktil, s5$mean_dist)))
# ---- Facebook ----
chk("FB reversal", f2(mean(fb$ktil)), f2(sd(fb$ktil)), f2(min(fb$ktil)), f2(max(fb$ktil)), f3(max(fb$ktil_se, na.rm = TRUE)), paste0("3 in ", sum(fb$krev == 3)), paste0("4 in ", sum(fb$krev == 4)))
chk("FB Q means", paste(f2(colMeans(fb[, paste0("Q", 1:6)])), collapse = "$, $"), f2(sd(fb$Q1)))
dif <- fb$ktil - fb$null_ktil; tt <- t.test(dif); chk("FB null", f3(mean(dif)), paste0(f3(tt$conf.int[1]), ",", f3(tt$conf.int[2])), f2(mean(fb$ktil)), f2(mean(fb$null_ktil)))
chk("FB MF", f3(cor(fb$ktil, fb$ktil_mf)), f2(sqrt(mean((fb$ktil_mf - fb$ktil)^2))), f2(mean(fb$ktil - fb$ktil_mf)), f3(cor(fb$ktil_mf, fb$null_ktil)), f2(mean(fb$ktil_mf - fb$null_ktil)))
fb$lnN <- log(fb$n); fb$mu10 <- fb$mean_deg / 10; fb$gap <- fb$ktil - fb$ktil_mf
for (mm in list(lm(ktil ~ cv + assort + clustering + mu10 + lnN, fb), lm(ktil ~ cv + assort + clustering + mu10 + lnN + modularity, fb), lm(gap ~ cv + assort + clustering + mu10 + lnN + modularity, fb))) { cf <- coef(summary(mm)); for (v in rownames(cf)) chk("Table 3", sprintf("%.3f", abs(cf[v, 1])), sprintf("(%.3f)", cf[v, 2])); chk("Table 3 R2", f3(summary(mm)$r.squared)) }
chk("FB shares", f2(mean(fb$c2)), f2(mean(fb$c3)), f2(mean(fb$c4)), f2(mean(fb$e5)), f2(mean(fb$e6)), f2(mean(fb$kP50 - fb$ktil)), f2(sd(fb$kP50 - fb$ktil)), f2(cor(fb$kP50, fb$ktil)))
# ---- Individual level ----
for (k in 1:5) { o <- ego[[paste0("own", k)]]; chk(paste("Table 4 P", k), f2(mean(o)), paste0("(", f2(sd(o)), ")"), paste0(f2(min(o)), "--", f2(max(o)))) }
chk("P3 networks", paste0("in ", sum(ego$own3 > 0.5), " of the 105"))
# ---- External networks ----
ex <- read.csv("../results/external_Z.csv"); ex$gap <- ex$ktil - ex$ktil_mf; ex$dn <- ex$ktil - ex$null_ktil
chk("external range", f2(min(ex$ktil)), f2(max(ex$ktil)), f3(max(ex$ktil_se, na.rm = TRUE)), paste0(min(ex$K), " to ", max(ex$K)))
chk("external corr", f2(cor(ex$null_ktil, ex$z)), f2(cor(ex$ktil_mf, ex$z)), f2(cor(ex$ktil, ex$z)), f3(cor(ex$ktil_mf, ex$null_ktil)), f2(max(abs(ex$ktil_mf - ex$null_ktil))))
chk("external offsets", f2(min(ex$dn)), f2(max(ex$dn)), f2(mean(ex$dn)))
chk("external modularity", f2(cor(ex$gap, ex$modularity, method = "spearman")), f2(cor(ex$gap, ex$modularity)), f2(cor(c(fb$ktil - fb$ktil_mf, ex$gap), c(fb$modularity, ex$modularity))), f2(cor(ex$gap, ex$clustering)), f2(abs(cor(ex$gap, ex$assort))))
chk("external kP50", f2(min(ex$kP50 - ex$ktil)), f2(max(ex$kP50 - ex$ktil)), f2(mean(ex$kP50 - ex$ktil)), f2(cor(ex$kP50, ex$ktil)))
for (i in seq_len(nrow(ex))) chk(paste("Table 4", ex$label[i]), paste0(f2(ex$ktil[i]), " & ", f2(ex$null_ktil[i]), " & ", f2(ex$ktil_mf[i])), paste0(f2(ex$modularity[i]), " & ", f2(ex$ktil[i])))
cat(sprintf("\naudit (with external): %d checks, %d missing\n", n, fails))

