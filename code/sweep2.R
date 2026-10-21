# Section 3 data: configuration-model sweep over n, mean degree, and degree CV.
src <- readLines("fb100_exact.R"); i <- grep("^res <- mclapply", src); eval(parse(text=src[1:(i-1)]))
NSRC <- 600; EXACT_N <- 1500
gen <- function(n, mu, cv, seed) { set.seed(seed)
  inv_r <- max(cv^2 - 1/mu, 1e-3); d <- rnbinom(n, size=1/inv_r, mu=mu); d[d<1] <- 1L; if (sum(d)%%2==1) d[1] <- d[1]+1L
  g <- simplify(sample_degseq(d, method="configuration")); cc <- components(g); induced_subgraph(g, which(cc$membership==which.max(cc$csize))) }
grid <- expand.grid(n=c(1000,4000,16000), mu=c(6,12,25,50), cv=c(0.6,0.94,1.4), rep=1:5)
res <- mclapply(seq_len(nrow(grid)), function(i) { p <- grid[i,]; tryCatch({
  g <- gen(p$n, p$mu, p$cv, seed=1000+i); deg <- degree(g); mu1 <- mean(deg); mu2 <- mean(deg^2)
  o <- shellstats(g, NSRC, 4, seed=i); r <- summ(o, mu1, mu2); mf <- mf_predict(deg, 4)
  nu <- mu2/mu1 - 1
  data.frame(p, n_lcc=vcount(g), mu1=mu1, cv_real=sd(deg)/mu1, nu=nu, x3=mu1*nu^2/vcount(g),
    t(setNames(r$B/mu1, paste0("Q",1:4))), t(setNames(mf/mu1, paste0("MFq",1:4))))
  }, error=function(e) NULL) }, mc.cores=min(4, detectCores()-1), mc.preschedule=FALSE)
write_csv(bind_rows(res), "sweep2.csv"); cat("DONE\n")
