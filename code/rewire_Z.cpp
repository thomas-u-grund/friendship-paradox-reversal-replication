// Degree-preserving double-edge swaps toward a target degree assortativity (Newman r) and/or global clustering (transitivity).
// Degrees are fixed, so the denominator of r and the number of connected triples are fixed; only sum(d_i d_j) over edges and the
// triangle count change, and both are updated incrementally for each proposed swap.
#include <Rcpp.h>
#include <unordered_set>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List rewire_targets(IntegerVector from, IntegerVector to, int n, double r_target, double t_target,
                    bool use_r, bool use_t, double w_r, double w_t, long max_prop, int seed, double temp = 0.0, IntegerVector block = IntegerVector::create()) {
  bool use_block = block.size() == n;
  int m = from.size(); std::vector<int> ea(m), eb(m), deg(n, 0);
  std::vector<std::unordered_set<int>> nb(n);
  for (int i = 0; i < m; i++) { ea[i] = from[i]; eb[i] = to[i]; deg[ea[i]]++; deg[eb[i]]++; nb[ea[i]].insert(eb[i]); nb[eb[i]].insert(ea[i]); }
  // fixed quantities
  double sx = 0, sx2 = 0, sxy = 0, triples = 0;
  for (int i = 0; i < m; i++) { sxy += (double)deg[ea[i]] * deg[eb[i]]; sx += 0.5 * (deg[ea[i]] + deg[eb[i]]); sx2 += 0.5 * ((double)deg[ea[i]] * deg[ea[i]] + (double)deg[eb[i]] * deg[eb[i]]); }
  for (int v = 0; v < n; v++) triples += 0.5 * deg[v] * (deg[v] - 1.0);
  double denom = sx2 / m - (sx / m) * (sx / m);
  auto assort = [&](double s) { return denom > 0 ? (s / m - (sx / m) * (sx / m)) / denom : 0.0; };
  // initial triangles
  double tri = 0;
  for (int i = 0; i < m; i++) { int a = ea[i], b = eb[i]; const auto &sa = nb[a].size() < nb[b].size() ? nb[a] : nb[b]; const auto &sb = nb[a].size() < nb[b].size() ? nb[b] : nb[a]; for (int x : sa) if (sb.count(x)) tri++; }
  tri /= 3.0;
  auto trans = [&](double t) { return triples > 0 ? 3.0 * t / triples : 0.0; };
  double r_cur = assort(sxy), t_cur = trans(tri), r0 = r_cur, t0 = t_cur;
  auto obj = [&](double r, double t) { double J = 0; if (use_r) J += w_r * (r - r_target) * (r - r_target); if (use_t) J += w_t * (t - t_target) * (t - t_target); return J; };
  double J_cur = obj(r_cur, t_cur);
  std::srand(seed); long acc = 0; Rcpp::RNGScope scope;
  auto common = [&](int u, int v, int ex1, int ex2) { // |N(u) ∩ N(v)| excluding vertices ex1, ex2
    const auto &su = nb[u].size() < nb[v].size() ? nb[u] : nb[v]; const auto &sv = nb[u].size() < nb[v].size() ? nb[v] : nb[u]; int c = 0;
    for (int x : su) if (x != ex1 && x != ex2 && sv.count(x)) c++; return c; };
  for (long it = 0; it < max_prop; it++) {
    int i = (int)(R::unif_rand() * m), j = (int)(R::unif_rand() * m); if (i == j) continue;
    int a = ea[i], b = eb[i], c = ea[j], d = eb[j]; if (R::unif_rand() < 0.5) std::swap(c, d);
    if (a == c || a == d || b == c || b == d) continue;            // four distinct vertices
    if (nb[a].count(d) || nb[c].count(b)) continue;                 // new edges must not exist
    if (use_block && !(block[a] == block[c] || block[b] == block[d])) continue;   // keeps the block-pair edge counts (planted modularity) fixed
    double dsxy = (double)deg[a] * deg[d] + (double)deg[c] * deg[b] - (double)deg[a] * deg[b] - (double)deg[c] * deg[d];
    double lost = common(a, b, -1, -1) + common(c, d, -1, -1);
    double gained = common(a, d, b, c) + common(c, b, a, d);       // in the graph without edges ab, cd
    double r_new = assort(sxy + dsxy), t_new = trans(tri - lost + gained), J_new = obj(r_new, t_new);
    bool ok = J_new < J_cur; if (!ok && temp > 0) ok = R::unif_rand() < std::exp(-(J_new - J_cur) / temp);
    if (!ok) continue;
    nb[a].erase(b); nb[b].erase(a); nb[c].erase(d); nb[d].erase(c); nb[a].insert(d); nb[d].insert(a); nb[c].insert(b); nb[b].insert(c);
    ea[i] = a; eb[i] = d; ea[j] = c; eb[j] = b; sxy += dsxy; tri += gained - lost; r_cur = r_new; t_cur = t_new; J_cur = J_new; acc++;
  }
  return List::create(_["from"] = IntegerVector(ea.begin(), ea.end()), _["to"] = IntegerVector(eb.begin(), eb.end()), _["r0"] = r0, _["t0"] = t0, _["r"] = r_cur, _["t"] = t_cur, _["accepted"] = (double)acc);
}
