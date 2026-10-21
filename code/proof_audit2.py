"""Second audit: behaviour of Qcal as lambda -> infinity, a sufficient condition for the crossing, and a denser grid for the single-crossing hypothesis.
Claim derived (see Appendix A): with m = H'(E*) < 1 and a(u) = E* - omega(u), a(nu u) = m a(u) (1 + o(1)), so
pi(d|s) = exp(-d omega(l s/(nu kappa))) - exp(-d omega(l s/kappa)) ~ d e^{-d E*} a(l s/kappa) (1/m - 1), and the s-dependent factor cancels in Qcal:
Q_inf = sum_d p_d d^2 e^{-d E*} / ( mu_1 sum_d p_d d e^{-d E*} ).
Run: python3 proof_audit2.py   (numpy, scipy); writes ../results/proof_audit2.txt"""
import numpy as np, sys
src = open("proof_audit.py").read(); exec(src.split("# ---------- T1")[0])
exec("def omega(L, u):\n    m = int(np.ceil(14 * np.log(10) / np.log(L['nu']))) + 5\n    E = np.asarray(u, float) / L['nu']**m\n    for _ in range(m): E = H(L, E)\n    return E")
exec(src[src.index("def Qcal"):src.index("def Qtilde")])
out = []
def log(s=""): print(s, flush=True); out.append(s)
def Qinf(L): Es = fixed_point(L); w = np.exp(-L["d"] * Es); return (L["p"] * L["d"]**2 * w).sum() / (L["mu1"] * (L["p"] * L["d"] * w).sum())
log("A. Limit of Qcal as lambda -> infinity: formula Q_inf vs Qcal at large lambda (largest lambda before 0/0 underflow), and crossing existence vs Q_inf < 1")
rows = []
for L in ALL:
    lam = np.exp(np.linspace(-8, 21, 59)); Q = np.array([Qcal(L, l) for l in lam]); ok = np.isfinite(Q); lamf, Qf = lam[ok], Q[ok]
    cross = int(np.sum(np.diff(np.sign(Qf - 1)) != 0)); qi = Qinf(L); rows.append((L["label"], L["nu"], cross, qi, Qf[-1], abs(Qf[-1] - qi), lamf[-1]))
dev = max(r[5] for r in rows); log(f"   70 laws: max |Qcal(largest finite lambda) - Q_inf| = {dev:.2e}")
agree = all((r[2] >= 1) == (r[3] < 1) for r in rows); log(f"   crossing exists  <=>  Q_inf < 1 in all 70 laws: {agree}  ({sum(r[3] < 1 for r in rows)} laws with Q_inf<1; {sum(r[3] >= 1 for r in rows)} with Q_inf>=1)")
log(f"   27 paper laws: Q_inf from {min(r[3] for r in rows[:27]):.3f} to {max(r[3] for r in rows[:27]):.3f}")
log("\nB. Denser grid (400 points of lambda in [e^-6, e^14]) for the 27 paper laws: number of crossings of Qcal=1, and whether Qcal is monotone decreasing")
bad = []
for L in laws:
    lam = np.exp(np.linspace(-6, 14, 400)); Q = np.array([Qcal(L, l) for l in lam]); ok = np.isfinite(Q); Q = Q[ok]
    c = int(np.sum(np.diff(np.sign(Q - 1)) != 0)); mono = bool(np.all(np.diff(Q) <= 1e-9))
    if c != 1 or not mono: bad.append((L["label"], c, mono))
log(f"   laws with exactly one crossing and Q monotone decreasing: {27 - len(bad)}/27" + (f"; exceptions {bad}" if bad else ""))
open("../results/proof_audit2.txt", "w").write("\n".join(out))
