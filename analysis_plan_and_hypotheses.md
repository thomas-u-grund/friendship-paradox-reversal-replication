# Analysis plan, pre-specified hypotheses and outcomes

The full plan (definitions, designs, decision rules, deviation log) follows `ANALYSIS_PLAN.md` below. The plan was written before the Stage 1 and Stage 2 data were analysed and is not externally registered; the file dates and the code are the only record. The community experiment (Stage 5) was designed after Stage 2; Stage 6, the coverage diagnostic and the proof audits were added after the main analysis.

## Hypotheses and outcomes (as in the paper)

- **H1.** Higher mean degree gives earlier reversal. *Supported* (-0.82 per log unit; Section 5.3).
- **H2.** Larger n gives later reversal. *Supported* (+0.36 per log unit; Section 5.3).
- **H3.** The reversal distance is  log n/ log nu+c(p), slope 1 within a degree law. *Supported with qualification:* fixed-law ratio mean 1.003 (range 0.87-1.16); the pooled test we specified gave 0.95, 0.97, and 1.23 for the three families, and the third family's design did not hold the degree law fixed (Section 5.4).
- **H4.** The effect of heterogeneity depends on density. *Supported:* earlier reversal at low density, later at high density, with the change-over density rising with n (Section 5.5).
- **H5.** Assortativity delays reversal. *Partly:* true at low clustering; with unrestricted rewiring it was reversed at clustering 0.20 and above (-0.16 and -0.15), but with the planted communities held fixed (Stage 6, added later) the effect is +0.14 and +0.21 at natural clustering and -0.03 and +0.03 at clustering 0.20 (Section 6.3).
- **H6.** Clustering delays reversal. *Supported,* including at fixed planted communities (Sections 6.2 and 6.3).
- **H7.** Real networks reverse later than their degree-preserving randomisations. *Supported* in 105 of 105 Facebook networks (Section 7.3) and in ten of ten public networks of other kinds (Section 7.5).
- **H8.** The individual-level crossing follows the reversal distance. *Supported* (lag 0.24-0.28 steps, correlation 0.98-0.99; Section 8).
- **H9.** Community structure delays reversal. *Supported* in rewiring (Section 6.3), in the Facebook networks (Section 7.4), and in the ten other networks (Section 7.5).

---

# Analysis plan for Version Z

*Fixed on 2026-10-20 before any new Z simulations. Changes after this date are recorded in the deviation log at the bottom, with a reason. The purpose is to separate hypotheses stated in advance from patterns found afterwards.*

## 1. Definitions

For a connected graph with n vertices, degrees d_v, mean degree μ₁ and second moment μ₂; shells N_k(v), sizes n_k(v); ν = μ₂/μ₁ − 1.

- **Pooled enrichment:** Q_k = B_k/μ₁, with B_k = Σ_v d_v n_k(v) / Σ_v n_k(v) (Theorem 1(a)); computed to the diameter.
- **Discrete reversal distance:** k_rev = min{k ≥ 2 : Q_k < 1}. Exists for non-regular graphs (Theorem 1(d)).
- **Continuous reversal distance (primary outcome):** k̃_rev = (k−1) + (Q_{k−1} − 1)/(Q_{k−1} − Q_k) at k = k_rev (linear interpolation in k). Robustness versions: interpolate log Q_k linearly; and the crossing of the cumulative curve Σ_{j≤k} E[n_j](Q_j − 1)/n (which is 0 at both ends by the conservation law).
- **Scaling variable:** z = log n / log ν. **Intercept:** c = k̃_rev − z.
- **Individual-level share:** P_k = share of people with n_k > 0 whose shell mean degree exceeds their own degree. **Persistence distance** k_P50 = the (interpolated) k at which P_k crosses 0.5.
- **Eligible-ego share** at k: share of vertices with n_k > 0. **Shell share:** E[n_k]/n.
- **Estimation:** exact (all vertices as sources) for n ≤ 4,000; otherwise 2,000 uniformly sampled sources with a bootstrap (200 resamples) SE for k̃_rev. All graphs restricted to the largest component. For sampled sources B_1 is set to μ₂/μ₁ exactly.

## 2. Hypotheses (directions stated in advance)

| ID | Hypothesis | Basis | Test |
|---|---|---|---|
| H1 | At fixed n and degree law, k̃_rev falls as mean degree rises | MF; B sweep | Stage 1, slope of k̃_rev on log(mean degree) at fixed n, CV |
| H2 | At fixed mean degree and CV, k̃_rev rises with log n | MF; B sweep | Stage 1 |
| H3 | k̃_rev = z + c(Λ) + error with slope 1 on z within each degree family | MF pilot: slope 0.996-1.006 | Stage 1: regress k̃_rev on z within family; slope 95% CI vs 1; report residual SD |
| H4 | Heterogeneity has two effects and its net sign depends on density: earlier crossing at low mean degree, later at high | MF pilot (mean degree 12: earlier; 50: later) | Stage 1: sign of ∂k̃_rev/∂CV at fixed (n, μ₁) in each density stratum; compare with MF |
| H5 | Higher degree assortativity delays reversal (later k̃_rev) at fixed degree sequence | B (residual +0.22 on Q₃) and FB regression (+0.65) | Stage 2 |
| H6 | Higher clustering delays reversal at fixed degree sequence | FB regression (+0.93, collinear) | Stage 2 (if clustering can be controlled with limited side effects) |
| H7 | Real networks reverse later than their degree-preserving randomisations | B: Q₃ 0.86 vs 0.70 | FB100: paired difference in k̃_rev, all 105 networks |
| H8 | The own-degree threshold below which individuals still have the paradox falls with k, and the 50% crossing of P_k occurs near k̃_rev | MF (one sign change per k verified on 5 networks); B terciles | MF + FB100 + Stage 1: correlation of k_P50 with k̃_rev; monotone threshold |
| H9 (secondary) | Stronger community structure (modularity) shifts k̃_rev beyond clustering and assortativity | proposal §6 | FB100 residual analysis; optional DC-SBM |

Direction in H9 is deliberately not specified; it is exploratory.

## 2b. Mean-field predictions recorded before the Stage 1 data are analysed

*(Written 2026-10-21, after the mean-field theory pilots `diagnostics/mf_theory_Z.R` and before any full Stage 1 result was looked at. Only a 20-graph feasibility pilot of Stage 1 had been seen. Stage 1 results will be compared with these predictions.)*

- **P1 (H3).** Within a degree family the continuous crossing distance is `z + c(Λ)` with slope 1 on z = log n / log ν. Model result: slope 0.989-1.003 for Negative-Binomial, lognormal and truncated power-law (γ = 2.5) degrees, n = 10³ to 10¹²; residual SD 0.003-0.026 (0.07-0.08 for the power law). Intercept c: Negative-Binomial 0.59-0.96, lognormal 0.54-0.89, power-law 0.75-1.02 (cases above), rising with CV for Negative-Binomial.
- **P2 (H4).** At n = 10⁴, raising CV from 0.4 to 2.0 changes the crossing by −1.78 (mean degree 6), −0.71 (12), −0.24 (25), −0.04 (50), +0.17 (100), +0.33 (200): earlier crossing at low mean degree, later at high, with the change of sign near mean degree 50 at this n. Finite-graph prediction: same signs; the sign-change density may move with n (to be checked).
- **P3 (H8).** In the model the 50% crossing of the individual-level share P_k (k_P50) lies 0.2-0.35 steps after the population crossing (mean degree 12: 3.39 vs 3.04, 4.72 vs 4.52, 6.26 vs 5.97, 7.72 vs 7.46 at n = 10³, 10⁵, 10⁷, 10⁹; mean degree 50: 2.46 vs 2.25, 3.47 vs 3.26, 4.47 vs 4.26, 5.47 vs 5.26), and both shift by the same amount per unit of log_ν n.
- **P4 (feasibility pilot, 20 graphs).** The mean-field crossing was within about 0.1 of the simulated mean in every pilot cell (for example 3.28 vs 3.31, 3.09 vs 3.16, 2.17 vs 2.20).

If Stage 1 contradicts any of P1-P4, report the contradiction; do not revise the predictions.

## 3. Designs

**Stage 1: degree sequence.**
- Degree law: Negative-Binomial with CV ∈ {0.4, 0.7, 1.0, 1.4}; sizes n ∈ {500, 1,000, 2,000, 4,000, 8,000, 16,000}; mean degree ∈ {6, 12, 25, 50, 100} (cells with mean degree > n/20 dropped). Configuration model (stub matching), simplified, largest component, minimum degree 1.
- ≥ 30 graphs per cell (≥ 60 for cells where k̃_rev's SD across graphs exceeds 0.1, decided from a first pass).
- Second and third degree families (lognormal; truncated power law with exponent 2.5) at n ∈ {2,000, 8,000}, mean degree ∈ {12, 50}, ≥ 30 graphs.
- Seeds: `set.seed(20000 + job index)`; job table stored with the outputs.
- Outputs per graph: n, mean degree, CV, ν, z, Q_1..Q_Δ, k_rev, k̃_rev (+ SE if sampled), P_1..P_Δ, k_P50, shell shares, eligible-ego shares, assortativity, transitivity.
- Analysis: within-family regressions for H1-H3; stratified comparison for H4; compare with the mean-field k̃_rev computed from each graph's realised degree sequence (agreement: RMSE, bias).

**Stage 2: beyond the degree sequence.**
- Base: Negative-Binomial sequences, n = 2,000, mean degree 25, CV 1.0 (and a second base: n = 4,000, mean degree 12, CV 1.4), each graph rewired from its own configuration-model draw so the degree sequence is exactly fixed.
- Assortativity targets {−0.3, −0.15, 0, +0.15, +0.3} by degree-preserving swaps accepted when they reduce |r − r_target| (as in B).
- Clustering targets: degree-preserving swaps accepted when they move the global clustering coefficient toward the target (starting from the unrewired value, up to the highest reachable value, recorded). Crossed with assortativity where both are reachable; otherwise one factor at a time.
- ≥ 30 graphs per cell; log, for every graph, the achieved assortativity, transitivity, modularity (Louvain), and mean path length. Feasibility is decided by a pilot (20 graphs per extreme target) before the full run: proceed only if the achieved range is ≥ 0.1 in transitivity and the side effects on the other factor are < 0.05.
- Analysis: regression of k̃_rev on achieved assortativity and transitivity with the side-effect covariates; report whether conclusions change if side effects are included.

**Facebook100.** 105 largest components, as in B; compute Q_k to the diameter, k_rev, k̃_rev with bootstrap SE (sampled networks), mean-field k̃_rev, degree-preserving null k̃_rev (one rewiring per network, 10 × edges swaps), modularity (Louvain), P_k and k_P50.

## 4. Statistical treatment

- Primary tests are H1-H3 (Stage 1) and H5, H6 (Stage 2); report effect sizes with 95% CIs, not only p-values. The 9 hypotheses are not corrected for multiplicity; the paper reports all of them, whether or not they hold.
- Facebook100 regressions use heteroskedasticity-robust SEs (HC3) and report variance inflation factors; they are labelled associations. Collinearity (clustering with log N, r = −0.87) is reported, and a model without clustering is shown alongside.
- Interpolation sensitivity (three crossing definitions) is reported in a table.
- Estimation noise: analyses of sampled networks use k̃_rev weighted by 1/SE² in a sensitivity analysis.

## 5. Decision rules (also in `ROADMAP.md`)

- **G1:** H3 supported (slope CI within [0.9, 1.1] in every degree family at n up to 16,000) → report the scaling statement; else report descriptive profiles only.
- **G2:** clustering/assortativity controllable (feasibility pilot passed) → report H5/H6 as controlled results; else report as correlates.
- **G3:** submit Z only if G1 passes and every number in the manuscript traces to a file in `data_out/` or `diagnostics/`.

## 6. Reporting commitments

Report all hypotheses (including those that fail); report the discrete k_rev alongside the continuous one; report eligible-ego shares and shell sizes at each k; give the seed table and code; state that the Facebook100 analysis is a validation on a single network class.

## 7. Deviation log

*2026-10-21. Stage 1 analysed (3,600 graphs; `an_stage1_Z.R`, `an_stage1_fe_Z.R`, `an_stage1_ratio_Z.R`).*

1. **Pre-specified H3/G1 test as written** (regress k̃_rev on z within degree family, slope 95% CI vs [0.9, 1.1]):
   Negative-Binomial 0.948 [0.940, 0.957] (pass), lognormal 0.973 [0.935, 1.011] (pass), truncated power law 1.225 [1.182, 1.268] (fail). By Negative-Binomial CV target: 0.989, 1.039, 1.074, 1.102 at CV 0.4, 0.7, 1.0, 1.4.
2. **Why the pre-specified test is not a clean test of the scaling statement (added after seeing the result).**
   - The statement is about a *fixed* degree law Λ with n varying; the pooled regression mixes cells with different Λ whose intercepts c(Λ) differ (Negative-Binomial: 0.38 to 1.05), and z uses each graph's realised ν, which varies between graphs of one cell.
   - The truncated-power-law cells let the degree cutoff grow with n (√(nμ)), so Λ is not fixed there; that family does not test the statement and its failure is not evidence against it.
3. **Refined tests (post hoc, reported alongside item 1):**
   - cell fixed effects on z: Negative-Binomial slope 1.043 [1.041, 1.045]; the mean-field values on the same graphs give 1.012;
   - **per-cell ratio** r = (d k̃_rev / d log n) × log ν within each (mean degree, CV) cell (20 Negative-Binomial cells, 4-6 sizes each): simulation mean 1.003 (SD 0.065, range 0.873-1.158; 17 of 20 cells within [0.9, 1.1]); mean-field on the same cells mean 0.961 (range 0.782-1.150); the two correlate 0.961 across cells.
   Interpretation: deviations from slope 1 at n = 10³-10⁴ are finite-size features that the mean-field model itself reproduces (larger than 1 in dense, low-CV cells; below 1 in dense, high-CV cells); the model's slope is 1.00 for n = 10³-10⁹ and for the other families to n = 10¹² (pilot). G1 is therefore judged **passed with qualification**: the scaling statement is a good first-order description, asymptotic in the model, with ±0.15 slope deviations at the sizes of real networks.
4. **Other Stage 1 results against predictions:** P4 confirmed (mean-field vs simulated crossing: corr 0.999, RMSE 0.093, mean-field later by 0.07, max |diff| 0.25); H1 and H2 confirmed (Negative-Binomial: −0.825 per log mean degree (−0.82 to two decimals), +0.355 per log n, R² 0.913 with CV fixed effects); H4 confirmed (change in crossing over CV 0.4 → 1.4: negative at mean degree 6 and 12 at every size, mixed at 25, positive at 50 and 100); P3 confirmed (k_P50 − k̃_rev mean 0.277, SD 0.148, correlation 0.984). Bootstrap SE of k̃_rev for sampled graphs: mean 0.007, max 0.015.

*2026-10-21 (evening). Stage 2 analysed (1,200 graphs; `an_stage2_Z.R`, `an_stage2_contrasts_Z.R`). Targets were reached: assortativity within 0.002 of target in every graph, clustering within 2.3e-4; degree sequences exactly fixed.*

5. **H6 (clustering delays reversal): supported as a total effect.** Continuous crossing, T = 0.30 vs the starting clustering (about 0.04), 95% CI: base A +0.340 [+0.327, +0.352] at assortativity 0 (+0.474 at −0.3; +0.206 at +0.3); base B +0.560 [+0.545, +0.575] (+0.710; +0.116). Regression coefficient on T: +1.38 (A) and +1.73 (B) per unit.
6. **H5 (assortativity delays reversal): only partly supported; the pre-stated direction fails at higher clustering.** Assortativity +0.3 vs −0.3: at the starting clustering +0.196 [+0.186, +0.206] (A) and +0.294 [+0.281, +0.308] (B); at T = 0.10 +0.086 and +0.023; at T = 0.20 **−0.158** and **−0.147**; at T = 0.30 −0.072 and **−0.300**. The assortativity × clustering interaction is −1.64 (A) and −3.56 (B). The effect of assortativity therefore depends on clustering, and H5 should be reported as conditional, not as a main effect.
7. **Side effects (G2).** Assortativity and clustering are independently controllable. But triangle-targeting swaps also raise **modularity** (T = 0.30 vs start: +0.344 [+0.331, +0.357] in A; +0.368 [+0.353, +0.383] in B) and **mean distance** (+0.164; +0.409). Correlation of clustering with modularity: 0.942 (A), 0.928 (B). Adding modularity to the regression changes clustering's coefficient from +1.38 to +0.20 (A) and from +1.73 to −1.01 (B); modularity (plus assortativity) alone explains R² 0.852 (A) and 0.890 (B). **The effect of clustering cannot be separated from community structure with this procedure; report it as clustering-with-community-structure.** G2 is passed with qualification.
8. **Consistency with real networks.** In cells with clustering 0.20 and assortativity 0 to 0.15 (close to the Facebook networks' averages), the crossing is later than the mean-field value by +0.07 to +0.17, similar in size to the offset observed in the 23 smallest Facebook networks (about +0.16). Suggestive only: that analysis is provisional until all 105 networks are done.

*2026-10-21 (night). WP4 analysed (105 Facebook100 largest components; `fb100_Z.R`, `an_fb100_Z.R`). A helper bug (`boot = 0` returned an empty list for sampled networks) stopped the first run after the 23 exact networks; fixed in `common_Z.R` (38 tests still pass) and rerun. Results below are for all 105.*

9. **Distribution of the reversal distance.** Discrete k_rev: 3 in 96, 4 in 8, 5 in 1. Continuous crossing: mean 2.752, SD 0.240, range 2.23-4.08; log-interpolation mean 2.708; bootstrap SE (82 sampled networks) mean 0.007, max 0.017. Eligible-ego share 1.00 through k = 4, 0.92 at k = 5, 0.49 at k = 6; mean shell share 0.33, 0.55, 0.10, 0.009 at k = 2..5.
10. **H7 supported.** Observed vs degree-preserving null (paired): 2.752 vs 2.549, difference +0.203 [0.194, 0.211]; observed later in 105/105 networks; correlation 0.985.
11. **Mean-field vs real networks.** Mean-field vs observed crossing: corr 0.982, RMSE 0.185, mean-field earlier by 0.179 (max 0.53). Mean-field vs the *null*: corr 0.999, mean-field later by only 0.024. So the mean-field model is accurate for a network with no structure beyond the degree sequence, and the real-minus-mean-field gap is real wiring.
12. **H3 across real networks** (cross-sectional; not the pre-specified fixed-law test): slope 1.425 [1.349, 1.502] on z, correlation 0.964, R² 0.929.
13. **Correlates of the crossing (OLS, N = 105).** Without modularity (R² 0.940): CV +0.172 (p = .016), assortativity +0.649 (p < .001), clustering +0.944 (p = .017), mean degree −0.063 per 10 (p < .001), log N +0.247 (p < .001). Adding modularity (R² 0.963): modularity +1.056 (p < .001); clustering falls to +0.394 (p = .21); assortativity +0.783; CV +0.257. Clustering is collinear with log N (r = −0.87) and only weakly with modularity here (r = −0.35).
14. **What explains the departure from mean-field** (ktil − ktil_mf; R² 0.705): assortativity +0.454 (p < .001), modularity +0.808 (p < .001), clustering +0.141 (n.s.), log N and CV n.s.
15. **H8 supported.** k_P50 mean 2.995; k_P50 − k̃_rev mean 0.243 (SD 0.059); correlation 0.990.
16. **Agreement across stages (descriptive).** Direction of the assortativity and modularity/clustering effects in the real networks matches Stage 2; the size of the real offset (+0.18 vs mean-field) is comparable to Stage 2 cells with clustering about 0.2. Not tested formally.

*2026-10-21 (night). WP5 and a WP1 result (decision by the author: separate clustering from communities, option 2; no extension). New scripts: `community_Z.R` (planted-community configuration model; block-preserving rewiring, a `block` argument added to `rewire_Z.cpp`), `stage5_Z.R`, `an_stage5_Z.R`; `mf_shift_check.R`.*

17. **WP5 design check (960 graphs; two bases, planted partitions of 10 and 20 blocks of 200 vertices; planted-community mixing 0.15/0.30/0.50/0.80, i.e. planted modularity 0.17-0.79; clustering targets start/0.10/0.20/0.30; 30 graphs per cell).** Block-preserving swaps keep the block-pair edge counts and degree sequence exactly: planted modularity changed by 0 in every graph; assortativity changed by at most 0.018 (at most 1e-15 in the strong-community cells). Clustering targets were reached except for the weakest communities (base A mix 0.8 target 0.30: reached 0.27; base B mix 0.8 targets 0.20 and 0.30: reached 0.15).
18. **Clustering has an effect at fixed planted community structure.** Continuous crossing, highest vs lowest clustering cell within each mixing level, 95% CI: base A +0.282 [0.272, 0.293] (planted Q 0.74), +0.194 [0.185, 0.203] (0.60), +0.224 [0.216, 0.231] (0.42), +0.253 [0.244, 0.262] (0.17); base B +0.639 [0.623, 0.655] (0.78), +0.455 [0.439, 0.471] (0.63), +0.373 [0.360, 0.386] (0.45), +0.175 [0.164, 0.186] (0.18). In the weak-community cells the detected (Louvain) modularity still rises with triangles (A: 0.17 → 0.50), so the separation is exact for the planted partition and approximate for detected communities; in the strong-community cells Louvain modularity barely moves (A mix 0.15: 0.74 → 0.75).
19. **Community strength has an effect at fixed clustering.** Strongest (mixing 0.15) minus weakest (0.80) at the same clustering target, 95% CI: base A +0.248 [0.238, 0.257] (T = 0.10), +0.265 [0.256, 0.275] (0.20), +0.306 [0.295, 0.317] (0.30; weak cell reached 0.27); base B +0.278 [0.267, 0.289] (0.10). (Base B at 0.20 and 0.30 are confounded: the weak cell reached only 0.15.) Regression of the crossing on achieved clustering and planted modularity: base A T +1.11 (se 0.03), Q +0.42 (se 0.012), R² 0.868; base B T +1.97 (se 0.05), Q +0.53 (se 0.02), R² 0.874; positive T × Q interaction (A +0.67, B +4.37).
20. **Mean distance accounts for almost all of the variation.** The crossing correlates 0.994 with the mean distance (0.991 within each base); adding mean distance to the regression gives R² 0.995 (A) and 0.992 (B) and reduces the clustering coefficient to +0.22 (A) and +0.31 (B) and the planted-modularity coefficient to about zero. Since the mean distance and the crossing are both summaries of how fast shells fill up, this is read as consistent with, not as an independent explanation of, the result: clustering and community structure delay reversal by slowing expansion.
21. **Net statement for the paper:** clustering and community strength each delay reversal, independently, by amounts of similar size (about +0.2 to +0.3 crossing steps for a change of 0.2 in clustering at fixed community structure, or of 0.57 in planted modularity at fixed clustering, in the bases studied), with a positive interaction. The Stage 2 finding that clustering and detected modularity cannot be separated in unrestricted rewiring stands; WP5 separates them by construction.
22. **WP1 result (mean-field translation property; `mf_shift_check.R`).** The mean-field recursion is autonomous in the cumulative exposure E_k and the frontier x_k: (E_k, x_{k+1}) = Φ(E_{k-1}, x_k), Φ(E, x) = (E + x, μ₁⁻¹ Σ_d p_d (d−1)(e^{-dE} − e^{-d(E+x)})), started at (0, s/(nμ₁)) for a source of degree s; n enters only through the starting value. Numerically, Q_{k+1}(νn) − Q_k(n) is 4e-4 at n = 10³, 1.6e-5 at 10⁵, 1.8e-7 at 10⁷, 6e-9 at 10⁹ (NB mean degree 12, CV 1.0; similar for mean degree 50 and for CV 1.4), over all shells holding between 1e-8 and 0.999 of the network. Origin of the result: Φ has a line of fixed points (E, 0), with Jacobian eigenvalues 1 and ν at the origin; the orbit leaves along the strong unstable manifold of the origin, so n shifts the orbit by log_ν n steps. Status: numerically verified; proof sketch in the manuscript appendix; a full proof (control of the neutral direction) not done.
23. **Monotone threshold (k = 1, 2).** m̃_k(s) is non-increasing in the ego degree s for k = 1 and k = 2 (likelihood-ratio ordering of the shell weights, because the cumulative exposure and the frontier are non-decreasing in s for k ≤ 2); hence m̃_k(s) − s has a single sign change and the paradox at distance k holds exactly for egos below a threshold degree. For k ≥ 3 the frontier need not be monotone in s; single crossing verified numerically on five networks only.

*2026-10-21 (late night). Additional work after the first manuscript draft.*

24. **Proof of Proposition 3 (mean-field translation property) completed.** Key observation: with H(E) = μ₁⁻¹ Σ_d p_d (d−1)(1 − e^{−dE}) the exposure recursion telescopes to the one-dimensional map E_{k+1} = δ_s + H(E_k), δ_s = s/(nμ₁) (verified against the two-variable recursion to 1e-15). With the Poincaré function ω of H at its repelling fixed point 0 (Koenigs; H(ω(u)) = ω(νu)), the orbit satisfies E_k = ω(ν^k δ_s/(ν−1)) + o(1) uniformly on compact scale ranges (Lemma A: a linear-phase estimate, using H(E) ≤ νE, plus the ν-Lipschitz property of H), so Q̃_k(n) = 𝒬(ν^k/n) + o(1). Consequences: translation by one step per factor ν in n; discrete reversal distance ⌈log_ν(nλ*)⌉; continuous reversal distance = u + π(u) + o(1) with an explicit continuous 1-periodic π. Numerical checks (`mf_proof_check.R`): λ* = 13.35 and c* = log_ν λ* = 0.828 for Negative-Binomial degrees, mean degree 12, CV 1.0; 𝒬(ν^k/n) reproduces Q̃_k(n) to 4e-4 (n = 10³), 1.7e-5 (10⁵), 1.9e-7 (10⁷); the interpolated crossing minus log_ν(nλ*) lies in [0, 0.016] and agrees with the π formula to 7e-4 at n = 10³ (differences vanish for larger n). Manuscript Appendix A now contains a full proof (cited theorem: Koenigs' linearisation, Milnor 2006, Thm 8.2: check).
25. **External networks (ten SNAP networks; `external_Z.R`, `an_external_Z.R`).** All ten reverse; observed continuous crossing 3.71 (Enron) to 11.36 (Amazon co-purchase), diameters seen 8-46, bootstrap SE ≤ 0.043. Degree-sequence benchmarks follow log n / log ν (corr 0.97 for the randomised networks, 0.98 for mean-field; observed 0.92). Mean-field matches the randomisations (corr 0.999; max difference 0.11). Observed later than both in 10/10, by 0.27 to 5.10 steps (mean 1.32). Offset (observed − mean-field) vs modularity: Spearman 0.99, Pearson 0.79 across the ten; Pearson 0.75 across all 115 networks; offset vs clustering 0.15, vs assortativity −0.06. k_P50 − k̃ = 0.74 to 4.84 (mean 1.69; corr 0.99). Caveat: ten networks from one repository, chosen for variety.
26. **Numerical audit (`audit_Z.R`).** 220 numbers quoted in the manuscript recomputed from the data files and checked against the text; one rounding error found and fixed (mean-degree coefficient −0.82, not −0.83); six formatting false alarms (negative numbers inside math delimiters, a share written as a percentage) checked by eye.
27. **Journal format.** `make_blinded.py` builds `main_Z_blinded.tex` with tables and figures after the references (30 pages).
