# Mechanism diagnostic: is enrichment Q_k a function of coverage R_k = sum_{j<=k} shell_share_j, whatever the wiring?
suppressMessages(library(dplyr)); suppressMessages(library(tidyr))
long <- function(f, cond) { d <- read.csv(f); d$id <- seq_len(nrow(d))
  Qs <- d %>% select(id, base, all_of(cond), clustering, assort, starts_with("Q")) %>% pivot_longer(starts_with("Q"), names_to = "k", values_to = "Q", names_prefix = "Q")
  Cs <- d %>% select(id, starts_with("c"), -clustering) %>% pivot_longer(starts_with("c"), names_to = "k", values_to = "c", names_prefix = "c")
  Cs <- Cs %>% mutate(k = as.integer(k)) %>% group_by(id) %>% arrange(k) %>% mutate(R = cumsum(ifelse(is.na(c), 0, c))) %>% ungroup()
  Qs %>% mutate(k = as.integer(k)) %>% inner_join(Cs, by = c("id", "k")) %>% filter(k >= 2, !is.na(Q), !is.na(c), R < 0.9995) }
r <- function(x) round(x, 3)
for (nm in c("stage5", "stage2")) {
  cond <- if (nm == "stage5") "mix" else "r_target"; f <- paste0("../results/", nm, "_full.csv")
  L <- long(f, cond); L$condition <- L[[cond]]
  cat("\n==", nm, ": shell rows k>=2 with R<1 :", nrow(L), "\n")
  res <- L %>% group_by(base) %>% group_modify(~{
    d <- .x
    m0 <- lm(Q ~ poly(R, 4), data = d)                       # Q as a flexible function of coverage only
    m1 <- lm(Q ~ poly(R, 4) + clustering + assort, data = d)   # + wiring
    mk <- lm(Q ~ factor(k), data = d); mkw <- lm(Q ~ factor(k) + clustering + assort, data = d)  # same but by shell index
    tibble(n = nrow(d), R2_shell_index = summary(mk)$r.squared, R2_index_plus_wiring = summary(mkw)$r.squared, R2_coverage = summary(m0)$r.squared, R2_coverage_plus_wiring = summary(m1)$r.squared,
           rmse_index = sqrt(mean(resid(mk)^2)), rmse_cov = sqrt(mean(resid(m0)^2)), rmse_cov_wiring = sqrt(mean(resid(m1)^2)),
           coef_clust_given_index = coef(mkw)["clustering"], coef_clust_given_cov = coef(m1)["clustering"],
           coef_assort_given_index = coef(mkw)["assort"], coef_assort_given_cov = coef(m1)["assort"]) }) %>% ungroup()
  print(res %>% mutate(across(where(is.numeric), r)), width = 250)
  write.csv(res, paste0("../results/coverage_", nm, ".csv"), row.names = FALSE) }
