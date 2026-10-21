# Stage 2 (ANALYSIS_PLAN.md): assortativity x clustering at fixed degree sequence. Usage: Rscript stage2_Z.R pilot|full [cores]
source("common_Z.R"); source("rewire_Z.R")
args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args) >= 1) args[1] else "pilot"; cores <- if (length(args) >= 2) as.integer(args[2]) else 12
outdir <- paste0("../results/stage2_", mode, "/"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
bases <- data.frame(base = c("A", "B"), n = c(2000, 4000), mu = c(25, 12), cv = c(1.0, 1.4), stringsAsFactors = FALSE)
r_t <- c(-0.30, -0.15, 0, 0.15, 0.30); t_t <- c(NA, 0.10, 0.20, 0.30)       # NA = hold clustering at its starting value
reps <- if (mode == "pilot") 2 else 30
grid <- expand.grid(base = bases$base, r_target = r_t, t_target = t_t, stringsAsFactors = FALSE)
if (mode == "pilot") grid <- grid[grid$base == "A" & grid$r_target %in% c(-0.3, 0.3) & grid$t_target %in% c(NA, 0.3), ]
jobs <- grid[rep(seq_len(nrow(grid)), each = reps), ]; jobs$rep <- rep(seq_len(reps), times = nrow(grid)); jobs$job <- seq_len(nrow(jobs)); jobs$seed <- 40000L + jobs$job
write_csv(jobs, paste0(outdir, "jobs.csv")); cat("mode", mode, ":", nrow(grid), "cells,", nrow(jobs), "graphs\n")
gen_base <- function(b, seed) { set.seed(seed); inv_r <- max(b$cv^2 - 1/b$mu, 1e-3); d <- rnbinom(b$n, size = 1/inv_r, mu = b$mu); d[d < 1] <- 1L; if (sum(d) %% 2 == 1) d[1] <- d[1] + 1L
  lcc(simplify(sample_degseq(d, method = "configuration"))) }
res <- mclapply(seq_len(nrow(jobs)), function(i) { j <- jobs[i, ]; f <- sprintf("%sjob_%05d.rds", outdir, j$job); if (file.exists(f)) return(NULL); t0 <- Sys.time()
  r <- tryCatch({ b <- bases[bases$base == j$base, ]; g0 <- gen_base(b, j$seed); r0 <- assortativity_degree(g0); t0g <- transitivity(g0)
    tt <- if (is.na(j$t_target)) t0g else j$t_target
    o <- rewire_to(g0, r_target = j$r_target, t_target = tt, max_prop = 8e6, w_r = 1, w_t = 1, seed = j$seed); g <- o$g
    cl <- cluster_louvain(g); rec <- graph_record(g, seed = j$seed, boot = 100)
    cbind(setNames(j, c("base", "r_target", "t_target", "rep", "job", "seed")), data.frame(r_start = r0, t_start = t0g, r_ach = o$info$r, t_ach = o$info$t, accepted = o$info$accepted, modularity = modularity(cl), n_comm = length(cl)), rec,
          mean_dist = sum(seq_len(10) * 0) + { p <- shell_profile(g); sum(seq_len(p$K) * colMeans(p$SZ)) / (vcount(g) - 1) }, secs = as.numeric(Sys.time() - t0, units = "secs")) }, error = function(e) { message("job ", j$job, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, f); NULL }, mc.cores = cores, mc.preschedule = FALSE)
fs <- list.files(outdir, pattern = "^job_.*\\.rds$", full.names = TRUE); out <- bind_rows(lapply(fs, readRDS)); write_csv(out, paste0("../results/stage2_", mode, ".csv")); cat("DONE", nrow(out), "of", nrow(jobs), "graphs\n")
