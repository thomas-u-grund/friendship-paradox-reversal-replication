# External networks (public SNAP edge lists; untrusted data: parsed only as two integer columns). Usage: Rscript external_Z.R [cores]
source("common_Z.R")
args <- commandArgs(trailingOnly = TRUE); cores <- if (length(args) >= 1) as.integer(args[1]) else 4
root <- "../external_networks/"; outdir <- "../results/external_Z/"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
nets <- data.frame(file = c("ca-GrQc", "ca-HepTh", "ca-CondMat", "ca-AstroPh", "email-Enron", "facebook_combined", "loc-brightkite_edges", "loc-gowalla_edges", "com-amazon", "com-dblp"),
  label = c("GrQc", "HepTh", "CondMat", "AstroPh", "Enron", "Facebook ego", "Brightkite", "Gowalla", "Amazon", "DBLP"),
  domain = c(rep("collaboration", 4), "email", "social (ego)", "location-based social", "location-based social", "co-purchase", "collaboration"), stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(nets)), function(i) { f <- paste0(root, "raw/", nets$file[i], ".txt.gz"); out <- paste0(outdir, nets$file[i], ".rds"); if (file.exists(out)) return(NULL); t0 <- Sys.time()
  r <- tryCatch({ e <- read.table(gzfile(f), comment.char = "#", colClasses = "numeric", col.names = c("a", "b")); stopifnot(ncol(e) == 2, all(is.finite(e$a)), all(is.finite(e$b)), all(e$a == round(e$a)))
    ids <- unique(c(e$a, e$b)); g <- lcc(simplify(graph_from_edgelist(cbind(match(e$a, ids), match(e$b, ids)), directed = FALSE))); n <- vcount(g)
    rec <- graph_record(g, seed = 21, boot = 200); cl <- cluster_louvain(g); set.seed(5); g0 <- lcc(rewire(g, keeping_degseq(niter = 10 * ecount(g)))); rec0 <- graph_record(g0, seed = 22, boot = 0, mf = FALSE)
    nul <- rec0[, c("ktil", "krev", "kP50", paste0("Q", 1:8))]; names(nul) <- paste0("null_", names(nul))
    cbind(data.frame(file = nets$file[i], label = nets$label[i], domain = nets$domain[i], m = ecount(g), max_deg = max(degree(g)), modularity = modularity(cl), n_comm = length(cl), secs = as.numeric(Sys.time() - t0, units = "secs")), rec, nul) },
    error = function(e) { message(nets$file[i], ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) saveRDS(r, out); NULL }, mc.cores = cores, mc.preschedule = FALSE)
fs <- list.files(outdir, pattern = "\\.rds$", full.names = TRUE); o <- bind_rows(lapply(fs, readRDS)); write_csv(o, "../results/external_Z.csv"); cat("DONE", nrow(o), "networks\n")
