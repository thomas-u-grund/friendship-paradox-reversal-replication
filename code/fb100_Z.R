# WP4: Facebook100 reversal distances. Usage: Rscript fb100_Z.R [cores]
source("common_Z.R")
args <- commandArgs(trailingOnly = TRUE); cores <- if (length(args) >= 1) as.integer(args[1]) else 4
root <- "../../"   # folder that contains facebook_networks/mtx_files/ (Facebook100 in Matrix Market format); outdir <- "../results/fb100_Z/"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(paste0(root, "facebook_networks/mtx_files/"), pattern = "\\.mtx$", full.names = TRUE); files <- files[file.info(files)$size < 15 * 1024 * 1024]; files <- files[order(-file.info(files)$size)]
res <- mclapply(files, function(f) { nm <- sub("\\.mtx$", "", basename(f)); out <- paste0(outdir, nm, ".rds"); if (file.exists(out)) return(NULL); t0 <- Sys.time()
  r <- tryCatch({ l <- readLines(f); dl <- l[!grepl("^%", l)][-1]; e <- read.table(text = dl)[, 1:2]; if (min(e) == 0) e <- e + 1
    g <- lcc(simplify(graph_from_edgelist(as.matrix(e), directed = FALSE))); rec <- graph_record(g, seed = 11, boot = 200); cl <- cluster_louvain(g)
    set.seed(3); g0 <- lcc(rewire(g, keeping_degseq(niter = 10 * ecount(g)))); rec0 <- graph_record(g0, seed = 12, boot = 0, mf = FALSE)
    nul <- rec0[, c("ktil", "ktil_log", "krev", "kP50", paste0("Q", 1:8), paste0("P", 1:8))]; names(nul) <- paste0("null_", names(nul))
    cbind(data.frame(network = nm, max_deg = max(degree(g)), modularity = modularity(cl), n_comm = length(cl), secs = as.numeric(Sys.time() - t0, units = "secs")), rec, nul) }, error = function(e) { message(nm, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, out); NULL }, mc.cores = cores, mc.preschedule = FALSE)
fs <- list.files(outdir, pattern = "\\.rds$", full.names = TRUE); o <- bind_rows(lapply(fs, readRDS)); write_csv(o, "../results/fb100_Z.csv"); cat("DONE", nrow(o), "networks\n")
