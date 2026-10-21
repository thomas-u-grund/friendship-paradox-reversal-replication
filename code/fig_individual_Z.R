# Figure 7: the paradox at the level of the individual (share P_k of people, by distance; and by own degree third). Input: ego_rel_all.csv (per-network shares).
r <- read.csv("../results/ego_rel_all.csv"); out <- "../figures/figure4_individual.pdf"
col <- c("#2a78d6","#eb6834","#1baf7a"); ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e3e2de"
pdf(out, width=7.2, height=3.3, family="Helvetica")
par(mfrow=c(1,2), mar=c(3.6,3.8,1.6,5.2), mgp=c(2.3,0.6,0), tcl=-0.25, cex.axis=0.8, cex.lab=0.9, col.axis=ink2, col.lab=ink, fg=ink2, las=1)
O <- as.matrix(r[,paste0("own",1:5)]); MF <- as.matrix(r[,paste0("mf_own",1:4)])
plot(NA, xlim=c(1,5), ylim=c(-0.04,1), xlab="Distance k", ylab="Share of people with paradox", xaxt="n", bty="l")
axis(1, at=1:5); abline(h=seq(0,1,0.25), col=grid, lwd=0.6); abline(h=0.5, col=ink2, lwd=1)
for (i in seq_len(nrow(O))) lines(1:5, O[i,], col=adjustcolor(col[1],0.18), lwd=0.7)
lines(1:5, colMeans(O), col=col[1], lwd=2.4); points(1:5, colMeans(O), pch=21, bg=col[1], col="white", cex=1.1)
lines(1:4, colMeans(MF), col=col[2], lwd=2.2, lty=2); points(1:4, colMeans(MF), pch=24, bg=col[2], col="white", cex=1.0)
legend("topright", legend=c("mean of networks (observed)", "mean-field prediction"), lty=c(1,2), lwd=c(2.4,2.2), col=c(col[1],col[2]), pch=c(21,24), pt.bg=c(col[1],col[2]), bty="o", box.col="white", bg="white", cex=0.68)
mtext("A", side=3, line=0.3, adj=-0.2, font=2, cex=1)
plot(NA, xlim=c(1,4), ylim=c(-0.04,1), xlab="Distance k", ylab="Share of people with paradox", xaxt="n", bty="l")
axis(1, at=1:4); abline(h=seq(0,1,0.25), col=grid, lwd=0.6); abline(h=0.5, col=ink2, lwd=1)
pch <- c(21,22,24)
for (tt in 1:3) { v <- sapply(1:4, function(k) mean(r[[paste0("own",k,"_t",tt)]])); lines(1:4, v, col=col[tt], lwd=2.2); points(1:4, v, pch=pch[tt], bg=col[tt], col="white", cex=1.1) }
text(4.2, c(0.50, 0.07, -0.02), c("lowest third", "middle third", "highest third"), adj=0, cex=0.68, col=ink, xpd=NA)
mtext("B", side=3, line=0.3, adj=-0.2, font=2, cex=1)
dev.off()
