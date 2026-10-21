suppressMessages(library(dplyr)); r <- read.csv("fb100_exact.csv"); Q <- as.matrix(r[,paste0("Q",1:6)])
krev <- apply(Q, 1, function(q) { k <- which(q[2:6] < 1)[1]; if (is.na(k)) NA else k+1 })
ktil <- apply(Q, 1, function(q) { k <- which(q[2:6] < 1)[1]; if (is.na(k)) NA else { kk <- k+1; kk - 1 + (q[kk-1]-1)/(q[kk-1]-q[kk]) } })
cat("discrete k_rev distribution (k>=2, from Q1..Q6):\n"); print(table(krev, useNA="ifany"))
cat("continuous (linearly interpolated) crossing: mean", round(mean(ktil,na.rm=T),3), "sd", round(sd(ktil,na.rm=T),3), "range", round(range(ktil,na.rm=T),2), "\n")
r$ktil <- ktil; r$lnN <- log(r$n); r$mu10 <- r$mu1/10
for (f in list(ktil ~ cv + assort + clustering + mu10 + lnN)) { m <- lm(f, r); cat("\nregression of continuous crossing distance: R2 =", round(summary(m)$r.squared,3), "\n"); print(round(coef(summary(m))[,1:4],3)) }
cat("\ncor(ktil, log n / log nu): nu = mu1*(1+cv^2)-1\n"); nu <- r$mu1*(1+r$cv^2)-1; z <- log(r$n)/log(nu); cat(round(cor(r$ktil, z, use="c"),3), " slope:", round(coef(lm(ktil~z,r))[2],3), "\n")
sw <- read.csv("sweep2.csv"); Qs <- as.matrix(sw[,paste0("Q",1:4)]); kk <- apply(Qs,1,function(q){k<-which(q[2:4]<1)[1]; if(is.na(k)) NA else k+1}); cat("\nconfig sweep (K=4 only): share with no crossing by k=4:", round(mean(is.na(kk)),3), "; table:\n"); print(table(kk,useNA="ifany"))
