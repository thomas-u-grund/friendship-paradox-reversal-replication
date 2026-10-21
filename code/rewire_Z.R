suppressMessages({library(Rcpp); library(igraph)}); sourceCpp("rewire_Z.cpp")
# g: igraph (simple, connected). Returns the rewired igraph and diagnostics. Vertex ids are 0-based inside C++.
rewire_to <- function(g, r_target = NA, t_target = NA, max_prop = 2e6, w_r = 1, w_t = 1, seed = 1, temp = 0, block = integer(0)) {
  E <- as_edgelist(g, names = FALSE) - 1L
  o <- rewire_targets(E[, 1], E[, 2], vcount(g), ifelse(is.na(r_target), 0, r_target), ifelse(is.na(t_target), 0, t_target), !is.na(r_target), !is.na(t_target), w_r, w_t, as.numeric(max_prop), seed, temp, as.integer(block))
  g2 <- make_graph(rbind(o$from, o$to) + 1L, n = vcount(g), directed = FALSE); list(g = g2, info = o[c("r0", "t0", "r", "t", "accepted")]) }
