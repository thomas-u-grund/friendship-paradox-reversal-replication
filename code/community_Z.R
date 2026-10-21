# WP5: planted-community configuration model and block-preserving rewiring.
source("common_Z.R"); source("rewire_Z.R")
# n vertices in B equal blocks; degrees from a given sequence; each vertex sends a fraction (1 - mix) of its stubs inside its block.
planted_graph <- function(deg, B = 10, mix = 0.3) {
  n <- length(deg); blk <- rep(seq_len(B), length.out = n); blk <- sample(blk)               # random assignment of vertices to blocks
  d_in <- rbinom(n, deg, 1 - mix); d_out <- deg - d_in
  edges <- list(); 
  for (b in seq_len(B)) { v <- which(blk == b); s <- rep(v, d_in[v]); if (length(s) %% 2 == 1) s <- s[-1]; s <- sample(s); if (length(s)) edges[[length(edges) + 1]] <- cbind(s[seq(1, length(s), 2)], s[seq(2, length(s), 2)]) }
  s <- rep(seq_len(n), d_out); if (length(s) %% 2 == 1) s <- s[-1]; s <- sample(s); if (length(s)) edges[[length(edges) + 1]] <- cbind(s[seq(1, length(s), 2)], s[seq(2, length(s), 2)])
  g <- simplify(make_graph(c(t(do.call(rbind, edges))), n = n, directed = FALSE)); keep <- which(components(g)$membership == which.max(components(g)$csize))
  list(g = induced_subgraph(g, keep), block = blk[keep]) }
planted_modularity <- function(g, block) modularity(g, membership = block)
