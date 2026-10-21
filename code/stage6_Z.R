# WP6 (added after external feedback): assortativity at fixed planted partition, fixed block-pair counts, fixed clustering. Usage: Rscript stage6_Z.R pilot|full [cores]
source("community_Z.R")
args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args) >= 1) args[1] else "pilot"; cores <- if (length(args) >= 2) as.integer(args[2]) else 12
outdir <- paste0("../results/stage6_", mode, "/"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
bases <- data.frame(base = c("A", "B"), n = c(2000, 4000), mu = c(25, 12), cv = c(1.0, 1.4), blocks = c(10, 20), stringsAsFactors = FALSE)
mixes <- c(0.8, 0.5, 0.15); t_t <- c(NA, 0.20); r_t <- c(-0.3, 0, 0.3); reps <- if (mode == "pilot") 2 else 30
grid <- expand.grid(base = bases$base, mix = mixes, t_target = t_t, r_target = r_t, stringsAsFactors = FALSE)
if (mode == "pilot") grid <- grid[grid$mix %in% c(0.5) & grid$t_target %in% c(NA, 0.2) & grid$r_target != 0, ]
jobs <- grid[rep(seq_len(nrow(grid)), each = reps), ]; jobs$rep <- rep(seq_len(reps), times = nrow(grid)); jobs$job <- seq_len(nrow(jobs)); jobs$seed <- 60000L + jobs$job
write_csv(jobs, paste0(outdir, "jobs.csv")); cat("mode", mode, ":", nrow(grid), "cells,", nrow(jobs), "graphs\n")
res <- mclapply(seq_len(nrow(jobs)), function(i) { j <- jobs[i, ]; f <- sprintf("%sjob_%05d.rds", outdir, j$job); if (file.exists(f)) return(NULL); t0 <- Sys.time()
  r <- tryCatch({ b <- bases[bases$base == j$base, ]; set.seed(j$seed); inv_r <- max(b$cv^2 - 1/b$mu, 1e-3); d <- rnbinom(b$n, size = 1/inv_r, mu = b$mu); d[d < 1] <- 1L
    p <- planted_graph(d, b$blocks, j$mix); g0 <- p$g; r0 <- assortativity_degree(g0); t0g <- transitivity(g0); Qp0 <- planted_modularity(g0, p$block); tt <- if (is.na(j$t_target)) t0g else j$t_target
    o <- rewire_to(g0, r_target = j$r_target, t_target = tt, max_prop = 8e6, seed = j$seed, block = p$block); g <- o$g; rec <- graph_record(g, seed = j$seed, boot = 100); pr <- shell_profile(g)
    cbind(setNames(j, c("base", "mix", "t_target", "r_target", "rep", "job", "seed")), data.frame(r_start = r0, r_ach = o$info$r, t_start = t0g, t_ach = o$info$t, planted_mod_start = Qp0, planted_mod = planted_modularity(g, p$block),
      louvain_mod = modularity(cluster_louvain(g)), mean_dist = sum(seq_len(pr$K) * colMeans(pr$SZ)) / (vcount(g) - 1), accepted = o$info$accepted), rec, secs = as.numeric(Sys.time() - t0, units = "secs")) },
    error = function(e) { message("job ", j$job, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, f); NULL }, mc.cores = cores, mc.preschedule = FALSE)
fs <- list.files(outdir, pattern = "^job_.*\\.rds$", full.names = TRUE); out <- bind_rows(lapply(fs, readRDS)); write_csv(out, paste0("../results/stage6_", mode, ".csv")); cat("DONE", nrow(out), "of", nrow(jobs), "graphs\n")
