suppressMessages(library(dplyr))
col <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"); ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e3e2de"
base_par <- function(...) par(mgp = c(2.2, 0.55, 0), tcl = -0.25, cex.axis = 0.78, cex.lab = 0.85, col.axis = ink2, col.lab = ink, fg = ink2, las = 1, ...)
se <- function(x) sd(x) / sqrt(length(x)); bars <- function(x, m, e, c) { ok <- is.finite(e) & e > 0; if (any(ok)) arrows(x[ok], (m - e)[ok], x[ok], (m + e)[ok], angle = 90, code = 3, length = 0.025, col = c, lwd = 0.9) }
s1 <- read.csv("../results/stage1_full.csv"); s2 <- read.csv("../results/stage2_full.csv"); s5 <- read.csv("../results/stage5_full.csv"); fb <- read.csv("../results/fb100_Z.csv")
# ---------- Figure 1: model vs simulation; scaling ----------
pdf("../figures/figureZ1_scaling.pdf", width = 7.2, height = 3.3, family = "Helvetica"); par(mfrow = c(1, 2), mar = c(3.6, 3.9, 1.6, 0.8)); base_par()
plot(NA, xlim = c(1.6, 6.2), ylim = c(1.6, 6.2), xlab = "Simulated crossing distance", ylab = "Mean-field crossing distance", bty = "l"); abline(h = 2:6, v = 2:6, col = grid, lwd = 0.6); abline(0, 1, col = ink2)
fam <- c(NB = col[1], LN = col[2], PL = col[3]); pchs <- c(NB = 21, LN = 24, PL = 22)
for (f in names(fam)) { d <- s1[s1$family == f, ]; points(d$ktil, d$ktil_mf, pch = pchs[f], bg = adjustcolor(fam[f], 0.35), col = adjustcolor(fam[f], 0.6), cex = 0.5, lwd = 0.3) }
legend("topleft", legend = c("Negative-Binomial", "lognormal", "power law"), pch = pchs, pt.bg = fam, col = "white", bty = "n", cex = 0.7); text(6.15, 1.75, sprintf("r = %.3f", cor(s1$ktil, s1$ktil_mf)), adj = 1, cex = 0.72, col = ink2)
mtext("A", side = 3, line = 0.3, adj = -0.2, font = 2, cex = 1)
nb <- s1[s1$family == "NB", ]; cells <- nb %>% group_by(mu_target, cv_target, n_target) %>% summarise(z = mean(z), k = mean(ktil), .groups = "drop")
plot(NA, xlim = range(cells$z) + c(-0.05, 0.1), ylim = range(cells$k) + c(-0.1, 0.1), xlab = expression("log n / log "*nu), ylab = "Simulated crossing distance", bty = "l"); abline(h = 2:6, v = seq(1, 4, 0.5), col = grid, lwd = 0.6)
cvs <- sort(unique(cells$cv_target)); for (i in seq_along(cvs)) for (m in sort(unique(cells$mu_target))) { d <- cells[cells$cv_target == cvs[i] & cells$mu_target == m, ] %>% arrange(n_target); if (nrow(d) < 2) next
  lines(d$z, d$k, col = adjustcolor(col[i], 0.8), lwd = 1.3); points(d$z, d$k, pch = 21, bg = col[i], col = "white", cex = 0.7); segments(d$z[1], d$k[1], d$z[nrow(d)], d$k[1] + (d$z[nrow(d)] - d$z[1]), col = adjustcolor(col[i], 0.5), lty = 3, lwd = 0.8) }
legend("topleft", legend = paste("CV", cvs), pch = 21, pt.bg = col[seq_along(cvs)], col = "white", lty = 1, lwd = 1.2, bty = "n", cex = 0.68, pt.cex = 0.8); text(max(cells$z), min(cells$k) - 0.02, "dotted: slope 1 from each line's first point", adj = 1, cex = 0.58, col = ink2)
mtext("B", side = 3, line = 0.3, adj = -0.2, font = 2, cex = 1); dev.off()
# ---------- Figure 2: heterogeneity effect changes sign with density ----------
pdf("../figures/figureZ2_heterogeneity.pdf", width = 5.2, height = 3.3, family = "Helvetica"); par(mar = c(3.6, 3.9, 1.4, 1.0)); base_par()
cm <- nb %>% group_by(n_target, mu_target, cv_target) %>% summarise(k = mean(ktil), .groups = "drop") %>% group_by(n_target, mu_target) %>% summarise(ch = k[cv_target == max(cv_target)] - k[cv_target == min(cv_target)], .groups = "drop")
plot(NA, xlim = c(5, 120), ylim = c(-1.35, 0.5), log = "x", xlab = "Mean degree (log scale)", ylab = "Change in crossing distance, CV 0.4 to 1.4", bty = "l"); abline(h = seq(-1, 0.5, 0.5), col = grid, lwd = 0.6); abline(h = 0, col = ink2)
for (i in 1:3) { n0 <- c(1000, 4000, 16000)[i]; d <- cm[cm$n_target == n0, ] %>% arrange(mu_target); lines(d$mu_target, d$ch, col = col[i], lwd = 1.8); points(d$mu_target, d$ch, pch = c(21, 22, 24)[i], bg = col[i], col = "white", cex = 1.1) }
legend("bottomright", legend = c("n = 1,000", "n = 4,000", "n = 16,000"), pch = c(21, 22, 24), pt.bg = col[1:3], col = "white", lty = 1, lwd = 1.4, bty = "n", cex = 0.72); dev.off()
# ---------- Figure 3: Stage 2 (assortativity x clustering) ----------
s2$tc <- ifelse(is.na(s2$t_target), "start", sprintf("%.2f", s2$t_target))
pdf("../figures/figureZ3_rewiring.pdf", width = 7.2, height = 3.2, family = "Helvetica"); par(mfrow = c(1, 3), mar = c(3.6, 3.9, 1.7, 0.6)); base_par()
for (b in c("A", "B")) { d <- s2[s2$base == b, ]; g <- d %>% group_by(r_target, tc) %>% summarise(T = mean(t_ach), k = mean(ktil), e = se(ktil), .groups = "drop")
  plot(NA, xlim = c(0, 0.33), ylim = range(g$k) + c(-0.06, 0.06), xlab = "Clustering coefficient", ylab = if (b == "A") "Crossing distance" else "", bty = "l"); abline(h = pretty(g$k), col = grid, lwd = 0.6)
  for (j in 1:3) { rt <- c(-0.3, 0, 0.3)[j]; h <- g[g$r_target == rt, ] %>% arrange(T); lines(h$T, h$k, col = col[j], lwd = 1.8); bars(h$T, h$k, h$e, col[j]); points(h$T, h$k, pch = c(21, 22, 24)[j], bg = col[j], col = "white", cex = 1) }
  if (b == "A") legend("topleft", legend = c("assort. -0.3", "assort. 0", "assort. +0.3"), pch = c(21, 22, 24), pt.bg = col[1:3], col = "white", lty = 1, lwd = 1.3, bty = "n", cex = 0.65, pt.cex = 0.85)
  mtext(sprintf("%s  Base %s", c(A = "A", B = "B")[b], b), side = 3, line = 0.3, adj = 0, font = 2, cex = 0.8) }
eff <- s2 %>% group_by(base, tc) %>% summarise(T = mean(t_ach), d = mean(ktil[r_target == 0.3]) - mean(ktil[r_target == -0.3]), e = sqrt(var(ktil[r_target == 0.3]) / sum(r_target == 0.3) + var(ktil[r_target == -0.3]) / sum(r_target == -0.3)), .groups = "drop")
plot(NA, xlim = c(0, 0.33), ylim = c(-0.4, 0.4), xlab = "Clustering coefficient", ylab = "Effect of assortativity (+0.3 minus -0.3)", bty = "l"); abline(h = seq(-0.4, 0.4, 0.2), col = grid, lwd = 0.6); abline(h = 0, col = ink2)
for (j in 1:2) { b <- c("A", "B")[j]; h <- eff[eff$base == b, ] %>% arrange(T); lines(h$T, h$d, col = col[j], lwd = 1.8); bars(h$T, h$d, 1.96 * h$e, col[j]); points(h$T, h$d, pch = c(21, 22)[j], bg = col[j], col = "white", cex = 1) }
legend("topright", legend = c("base A", "base B"), pch = c(21, 22), pt.bg = col[1:2], col = "white", lty = 1, lwd = 1.3, bty = "n", cex = 0.68); mtext("C  Assortativity effect", side = 3, line = 0.3, adj = 0, font = 2, cex = 0.8); dev.off()
# ---------- Figure 4: WP5 (clustering at fixed planted communities) ----------
s5$tc <- ifelse(is.na(s5$t_target), "start", sprintf("%.2f", s5$t_target))
pdf("../figures/figureZ4_communities.pdf", width = 7.2, height = 3.2, family = "Helvetica"); par(mfrow = c(1, 2), mar = c(3.6, 3.9, 1.7, 0.8)); base_par()
for (b in c("A", "B")) { d <- s5[s5$base == b, ]; g <- d %>% group_by(mix, tc) %>% summarise(T = mean(t_ach), Q = mean(planted_mod), k = mean(ktil), e = se(ktil), .groups = "drop")
  plot(NA, xlim = c(0, 0.33), ylim = range(g$k) + c(-0.08, 0.08), xlab = "Clustering coefficient", ylab = if (b == "A") "Crossing distance" else "", bty = "l"); abline(h = pretty(g$k), col = grid, lwd = 0.6)
  mixes <- sort(unique(g$mix)); for (j in seq_along(mixes)) { h <- g[g$mix == mixes[j], ] %>% arrange(T); lines(h$T, h$k, col = col[j], lwd = 1.8); bars(h$T, h$k, h$e, col[j]); points(h$T, h$k, pch = c(21, 22, 24, 23)[j], bg = col[j], col = "white", cex = 1)
    text(max(h$T) + 0.005, h$k[nrow(h)], sprintf("Q = %.2f", mean(h$Q)), adj = 0, cex = 0.6, col = ink, xpd = NA) }
  mtext(sprintf("%s  Base %s (planted communities)", b, b), side = 3, line = 0.3, adj = 0, font = 2, cex = 0.8) }
dev.off()
# ---------- Figure 5: Facebook100 ----------
pdf("../figures/figureZ5_facebook.pdf", width = 7.2, height = 3.3, family = "Helvetica"); par(mfrow = c(1, 2), mar = c(3.6, 3.9, 1.6, 0.8)); base_par()
rg <- range(c(fb$ktil, fb$ktil_mf, fb$null_ktil)); plot(NA, xlim = rg, ylim = rg, xlab = "Observed crossing distance", ylab = "Crossing distance of the comparison model", bty = "l"); abline(h = pretty(rg), v = pretty(rg), col = grid, lwd = 0.6); abline(0, 1, col = ink2)
points(fb$ktil, fb$null_ktil, pch = 24, bg = adjustcolor(col[2], 0.75), col = "white", cex = 0.9); points(fb$ktil, fb$ktil_mf, pch = 1, col = col[1], cex = 1.15, lwd = 1.1)
legend("topleft", legend = c("mean-field model", "degree-preserving rewiring"), pch = c(1, 24), pt.bg = c(NA, col[2]), col = c(col[1], "white"), bty = "n", cex = 0.7, pt.lwd = 1.1); mtext("A", side = 3, line = 0.3, adj = -0.2, font = 2, cex = 1)
fb$gap <- fb$ktil - fb$ktil_mf; m <- lm(gap ~ assort + modularity, fb); rg2 <- range(c(fb$gap, fitted(m)))
plot(fitted(m), fb$gap, xlab = "Fitted from assortativity and modularity", ylab = "Observed minus mean-field crossing", xlim = rg2, ylim = rg2, pch = 21, bg = adjustcolor(col[1], 0.7), col = "white", cex = 0.9, bty = "l"); abline(h = pretty(rg2), v = pretty(rg2), col = grid, lwd = 0.6); abline(0, 1, col = ink2)
points(fitted(m), fb$gap, pch = 21, bg = adjustcolor(col[1], 0.7), col = "white", cex = 0.9); text(rg2[2], rg2[1] + 0.03, sprintf("R2 = %.2f", summary(m)$r.squared), adj = 1, cex = 0.72, col = ink2); mtext("B", side = 3, line = 0.3, adj = -0.2, font = 2, cex = 1); dev.off()
cat("figures written; Facebook gap model R2:", round(summary(m)$r.squared, 3), "\n")
# ---------- Figure 6: external networks ----------
ex <- read.csv("../results/external_Z.csv"); ex$gap <- ex$ktil - ex$ktil_mf; fb$gap <- fb$ktil - fb$ktil_mf
pdf("../figures/figureZ6_external.pdf", width = 7.2, height = 3.4, family = "Helvetica"); par(mfrow = c(1, 2), mar = c(3.6, 3.9, 1.6, 0.8)); base_par()
rg <- c(2, 12); plot(NA, xlim = rg, ylim = rg, xlab = "Observed crossing distance", ylab = "Degree-preserving randomisation", bty = "l"); abline(h = seq(2, 12, 2), v = seq(2, 12, 2), col = grid, lwd = 0.6); abline(0, 1, col = ink2)
points(fb$ktil, fb$null_ktil, pch = 21, bg = adjustcolor(ink2, 0.35), col = "white", cex = 0.7); points(ex$ktil, ex$null_ktil, pch = 21, bg = col[1], col = "white", cex = 1.2)
lab <- ex[ex$gap > 0.8 | ex$label %in% c("Enron"), ]; text(lab$ktil, lab$null_ktil, lab$label, pos = 1, cex = 0.6, col = ink, offset = 0.45)
legend("topleft", legend = c("Facebook100 (105)", "other networks (10)"), pch = 21, pt.bg = c(adjustcolor(ink2, 0.35), col[1]), col = "white", bty = "n", cex = 0.7); mtext("A", side = 3, line = 0.3, adj = -0.2, font = 2, cex = 1)
plot(NA, xlim = c(0.4, 1), ylim = c(0, 5.4), xlab = "Modularity (Louvain)", ylab = "Observed minus mean-field crossing", bty = "l"); abline(h = 0:5, col = grid, lwd = 0.6); abline(h = 0, col = ink2)
points(fb$modularity, fb$gap, pch = 21, bg = adjustcolor(ink2, 0.35), col = "white", cex = 0.7); points(ex$modularity, ex$gap, pch = 21, bg = col[1], col = "white", cex = 1.2)
sel <- ex$label %in% c("Amazon", "GrQc", "Facebook ego", "DBLP", "HepTh", "Enron"); text(ex$modularity[sel], ex$gap[sel], ex$label[sel], pos = ifelse(ex$label[sel] %in% c("GrQc", "HepTh", "Enron"), 2, 4), cex = 0.6, col = ink, offset = 0.4, xpd = NA); mtext("B", side = 3, line = 0.3, adj = -0.2, font = 2, cex = 1); dev.off(); cat("figure 6 written\n")
