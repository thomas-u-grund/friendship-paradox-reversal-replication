# Replication pack: "When Does the Friendship Paradox Reverse? Network Structure and the Limits of Degree Bias"

This folder contains the code and the result files behind every number, table and figure in the paper. It does not contain the Facebook100 or SNAP network data (see below).

## Layout
- `code/`: R and C++ scripts (run from this folder; paths are relative, `../results/`, `../figures/`, and the manuscript at `../manuscript.tex`).
- `results/`: result files produced by the scripts, one row per simulated or empirical network. Everything in the paper's tables can be recomputed from these without re-running the simulations.
- `figures/`: the figure PDFs used in the paper.
- `external_networks/get_external_networks.sh`: downloads the ten public SNAP edge lists.

## Software
R (4.x) with `igraph`, `dplyr`, `tidyr`, `readr`, `Rcpp`, `parallel`; a C++ compiler for `rewire_Z.cpp` (compiled on the fly by `Rcpp::sourceCpp`); Python 3 with `numpy` and `scipy` for `proof_audit.py` only. Simulations used up to 12 cores; seeds are fixed in the job tables, so reruns reproduce the results.

## Data not included
- **Facebook100** (Traud, Mucha, & Porter, 2012): obtain from the original distribution and place the `.mat` files in a folder named `facebook_networks/` two levels above `code/` (the path is set at the top of `fb100_Z.R`). Not redistributed here.
- **SNAP networks** (Leskovec & Krevl, 2014): `sh external_networks/get_external_networks.sh`. The script `external_Z.R` expects them in `raw/`.

## From script to result
| Script | Produces | Used for |
|---|---|---|
| `common_Z.R` (tests: `test_common_Z.R`) | shell profiles, reversal distances, mean-field model | all |
| `stage1_Z.R`, `an_stage1_Z.R`, `an_stage1_fe_Z.R`, `an_stage1_ratio_Z.R` | `results/stage1_full.csv` (3,600 configuration-model graphs) | Section 5, Table 1, Figures 2 and 3 |
| `stage2_Z.R` with `rewire_Z.cpp`, `an_stage2_Z.R`, `an_stage2_contrasts_Z.R` | `results/stage2_full.csv` (1,200 graphs) | Section 6.2, Table 2, Figure 4 |
| `stage5_Z.R`, `community_Z.R`, `an_stage5_Z.R` | `results/stage5_full.csv` (960 graphs) | Section 6.3, Table 2, Figure 5 |
| `stage6_Z.R`, `an_stage6_Z.R`, `an_assort_mod.R` | `results/stage6_full.csv`, `stage6_effects.csv` (1,080 graphs; added after the main analysis) | Section 6.3, Table 2 |
| `an_coverage_Z.R`, `an_coverage2_Z.R` | `results/coverage_*.csv` | Section 6.3 (share of the network reached at the crossing) |
| `fb100_Z.R`, `an_fb100_Z.R`, `an_ego.R` | `results/fb100_Z.csv` (105 networks) | Sections 7-8, Tables 3 and 5, Figures 6 and 8 |
| `external_Z.R`, `an_external_Z.R` | `results/external_Z.csv` (10 networks) | Section 7.5, Table 4, Figure 7 |
| `an_review_checks.R` | clustering regime of the Facebook networks, leave-one-out for the ten networks | Sections 6.3, 7.4, 7.5 |
| `mf_proof_check.R`, `mf_single_crossing.R`, `mf_threshold_check.R`, `proof_audit.py`, `proof_audit2.py` | `results/mf_*.csv`, `proof_audit.txt` | Appendices A and B (numerical checks of the proofs) |
| `fig_schematic_Z.R` | `figures/figure_intro_schematic.pdf` | Figure 1 (illustrative network) |
| `fig_individual_Z.R` | `figures/figure4_individual.pdf` | Figure 8 |
| `tables_Z.R`, `fig_Z.R` | tables and figures | Figures 2-7 and the tables |
| `audit_Z.R` | recomputes the numbers quoted in the manuscript and checks that each appears in the text | consistency check |

## Pre-specification and departures
The analysis plan, the hypotheses H1-H9 and the mean-field predictions were written before the Stage 1 data were analysed (see `analysis_plan_and_hypotheses.md`, which has the hypotheses, outcomes and departures). Four analyses were added later and are labelled as such in the paper: the pooled slope test for the power-law family departs from its pre-specified form, the coverage-at-the-crossing diagnostic, the assortativity experiment inside planted partitions (Stage 6), and the proof audits.

## Checking the mathematics
`proof_audit.py` tests each claim of Appendix A numerically for 27 degree laws used in the paper and 43 random degree laws, including near-critical ones: properties of H, the Poincare function, Lemma A (error against delta), the limit profile of Q, the translation error, the discrete and continuous reversal-distance statements with n close to integer crossings, and the existence of the crossing. Output: `results/proof_audit.txt`.

## Use of AI tools
Claude Sonnet 5.5 (Anthropic, October 2026) assisted with the programming and with preparing these replication materials; the author reviewed the output and takes responsibility for it.
