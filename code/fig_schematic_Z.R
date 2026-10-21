# Illustrative figure (Figure 1): a small network drawn in rings by distance from one person, and the enrichment Q_k of each distance shell.
suppressMessages(library(igraph))
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e3e2de"; blue <- "#2a78d6"; orange <- "#eb6834"
seed <- 16; set.seed(seed); n <- 36; mu <- 3.0; cv <- 1.2; inv_r <- max(cv^2 - 1/mu, 1e-3)
d <- rnbinom(n, size = 1/inv_r, mu = mu); d[d < 1] <- 1L; if (sum(d) %% 2) d[1] <- d[1] + 1L
g <- simplify(sample_degseq(d, method = "vl")); stopifnot(is_connected(g))
deg <- degree(g); D <- distances(g); mu1 <- mean(deg); K <- max(D)
Q <- sapply(1:K, function(k) { M <- (D == k); sum(M %*% deg) / sum(M) / mu1 })
# focal person: a vertex of degree 3 whose shells reach the full range
cand <- which(deg == 3); focal <- cand[which.max(sapply(cand, function(v) max(D[v, ])))]
dist_f <- D[focal, ]; nk <- table(dist_f)
# ring layout: shell k on radius k, nodes ordered by the mean angle of their parents
ang <- numeric(n); ang[focal] <- 0; xy <- matrix(0, n, 2)
sh1 <- which(dist_f == 1); ang[sh1] <- 2 * pi * (seq_along(sh1) - 0.5) / length(sh1)
for (k in 2:max(dist_f)) { v <- which(dist_f == k); pa <- sapply(v, function(x) { p <- intersect(neighbors(g, x), which(dist_f == k - 1)); atan2(mean(sin(ang[p])), mean(cos(ang[p]))) })
  o <- order(pa %% (2 * pi)); ang[v[o]] <- 2 * pi * (seq_along(v) - 0.5) / length(v) + (k %% 2) * 0.15 }
xy[, 1] <- dist_f * cos(ang); xy[, 2] <- dist_f * sin(ang); xy[focal, ] <- 0
nodecol <- ifelse(dist_f >= 1 & Q[pmax(dist_f, 1)] >= 1, blue, orange)   # same colours as the bars: shell above average = blue, below = orange
pdf("../figures/figure_intro_schematic.pdf", width = 7.2, height = 3.5, family = "Helvetica")
layout(matrix(1:2, 1), widths = c(1.1, 1)); par(mar = c(1.2, 0.6, 2.0, 0.6), mgp = c(2.3, 0.6, 0))
plot(NA, xlim = range(xy[, 1]) * 1.08, ylim = c(min(xy[, 2]) * 1.08, max(xy[, 2]) * 1.16), asp = 1, axes = FALSE, xlab = "", ylab = "")
for (k in 1:max(dist_f)) symbols(0, 0, circles = k, inches = FALSE, add = TRUE, fg = grid, lwd = 0.8)
E <- as_edgelist(g, names = FALSE); segments(xy[E[, 1], 1], xy[E[, 1], 2], xy[E[, 2], 1], xy[E[, 2], 2], col = adjustcolor(ink2, 0.45), lwd = 0.8)
cx <- 0.45 + 0.30 * sqrt(deg); points(xy[, 1], xy[, 2], pch = 21, bg = nodecol, col = "white", cex = cx, lwd = 0.8)
points(xy[focal, 1], xy[focal, 2], pch = 21, bg = ink, col = "white", cex = cx[focal] + 0.5, lwd = 0.8)
text(0, -0.5, "you", cex = 0.62, col = ink, font = 2)
for (k in 1:max(dist_f)) text(0.06 * k, k, bquote(k == .(k)), cex = 0.6, col = ink2, pos = 3, offset = 0.2)
mtext("A  One network, by distance from one person", side = 3, line = 0.2, adj = 0, font = 2, cex = 0.8)
mtext("circle size = number of friends; blue = better connected than average, orange = less", side = 1, line = 0.1, adj = 0.5, cex = 0.58, col = ink2)
par(mar = c(3.9, 4.6, 2.0, 0.8), mgp = c(2.4, 0.6, 0), tcl = -0.25, cex.axis = 0.8, cex.lab = 0.85, col.axis = ink2, col.lab = ink, fg = ink2, las = 1)
bp <- barplot(Q, names.arg = 1:K, col = ifelse(Q >= 1, blue, orange), border = NA, ylim = c(0, max(Q) * 1.15), xlab = "Distance k from a person", ylab = "")
mtext("Mean degree at distance k,\nrelative to everyone (Q_k)", side = 2, line = 2.6, las = 0, cex = 0.72, col = ink)
abline(h = 1, col = ink, lty = 2, lwd = 1)
text(bp[K] + 0.4, 1.07, "everyone (= 1)", adj = 1, cex = 0.62, col = ink)
kr <- which(Q < 1)[1]; text(bp[kr] + 0.45, 1.13, "reversal", cex = 0.7, col = orange, font = 2, adj = 0.5)
text(bp[1], Q[1] + 0.07, sprintf("%.2f", Q[1]), cex = 0.6, col = ink); text(bp[2], Q[2] + 0.07, sprintf("%.2f", Q[2]), cex = 0.6, col = ink)
mtext("B  Mean degree at each distance", side = 3, line = 0.2, adj = 0, font = 2, cex = 0.8)
dev.off()
cat("Q_k:", round(Q, 2), "| n", n, "edges", ecount(g), "| focal degree", deg[focal], "| shell sizes from focal:", paste(nk[-1], collapse = " "), "\n")
