# Stage 1 (ANALYSIS_PLAN.md): degree sequence and reversal distance. Usage: Rscript stage1_Z.R pilot|full [cores]
source("common_Z.R")
args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args) >= 1) args[1] else "pilot"; cores <- if (length(args) >= 2) as.integer(args[2]) else 12
outdir <- paste0("../results/stage1_", mode, "/"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- degree-sequence generators (minimum degree 1, even sum) ----
fix_seq <- function(d) { d <- pmax(as.integer(d), 1L); if (sum(d) %% 2 == 1) d[1] <- d[1] + 1L; d }
gen_degrees <- function(family, n, mu, cv) switch(family,
  NB = { inv_r <- max(cv^2 - 1/mu, 1e-3); fix_seq(rnbinom(n, size = 1/inv_r, mu = mu)) },
  LN = { s2 <- log(1 + cv^2); fix_seq(round(rlnorm(n, log(mu) - s2/2, sqrt(s2)))) },
  PL = { gam <- 2.5; dmax <- floor(sqrt(n * mu)); best <- NULL
         for (dmin in 1:max(1, floor(mu))) { d <- dmin:dmax; p <- d^-gam; p <- p/sum(p); m <- sum(d*p); if (is.null(best) || abs(m - mu) < best$e) best <- list(e = abs(m - mu), d = d, p = p) }
         fix_seq(sample(best$d, n, replace = TRUE, prob = best$p)) })
make_graph <- function(family, n, mu, cv) lcc(simplify(sample_degseq(gen_degrees(family, n, mu, cv), method = "configuration")))

# ---- job table ----
nb <- expand.grid(family = "NB", n = c(500, 1000, 2000, 4000, 8000, 16000), mu = c(6, 12, 25, 50, 100), cv = c(0.4, 0.7, 1.0, 1.4), stringsAsFactors = FALSE)
nb <- nb[nb$mu <= nb$n / 20, ]
ln <- expand.grid(family = "LN", n = c(2000, 8000), mu = c(12, 50), cv = c(0.7, 1.4), stringsAsFactors = FALSE)
pl <- expand.grid(family = "PL", n = c(2000, 8000), mu = c(12, 50), cv = NA_real_, stringsAsFactors = FALSE)
cells <- rbind(nb, ln, pl)
reps <- if (mode == "pilot") 2 else 30
if (mode == "pilot") cells <- cells[(cells$family == "NB" & cells$n %in% c(1000, 4000) & cells$mu %in% c(12, 50) & cells$cv %in% c(0.7, 1.4)) | (cells$family == "LN" & cells$n == 2000 & cells$mu == 12 & cells$cv == 0.7) | (cells$family == "PL" & cells$n == 2000 & cells$mu == 12), ]
jobs <- cells[rep(seq_len(nrow(cells)), each = reps), ]; jobs$rep <- rep(seq_len(reps), times = nrow(cells)); jobs$job <- seq_len(nrow(jobs)); jobs$seed <- 20000L + jobs$job
write_csv(jobs, paste0(outdir, "jobs.csv")); cat("mode", mode, ":", nrow(cells), "cells,", nrow(jobs), "graphs\n")

# big jobs first for load balance; each job saved separately so a killed run can resume
ord <- order(-jobs$n * jobs$mu)
res <- mclapply(ord, function(i) { j <- jobs[i, ]; f <- sprintf("%sjob_%05d.rds", outdir, j$job); if (file.exists(f)) return(NULL)
  t0 <- Sys.time(); set.seed(j$seed)
  r <- tryCatch({ g <- make_graph(j$family, j$n, j$mu, j$cv); cbind(setNames(j, c("family","n_target","mu_target","cv_target","rep","job","seed")), graph_record(g, seed = j$seed), secs = as.numeric(Sys.time() - t0, units = "secs")) }, error = function(e) NULL)
  if (!is.null(r)) saveRDS(r, f); NULL }, mc.cores = cores, mc.preschedule = FALSE)
fs <- list.files(outdir, pattern = "^job_.*\\.rds$", full.names = TRUE); out <- bind_rows(lapply(fs, readRDS))
write_csv(out, paste0("../results/stage1_", mode, ".csv")); cat("DONE", nrow(out), "of", nrow(jobs), "graphs\n")
