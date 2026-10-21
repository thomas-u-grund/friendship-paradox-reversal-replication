"""Numerical audit of Appendix A (proof of Proposition 3), claim by claim.
Run: python3 proof_audit.py   (numpy, scipy). Writes ../results/proof_audit.txt.
Each test states the claim from the manuscript, then reports the worst case found over a set of degree laws:
27 'paper' laws (NB, lognormal, truncated power law) plus random hostile laws."""
import numpy as np, sys
from scipy import stats, optimize
rng = np.random.default_rng(20261021)
out = []
def log(s=""): print(s); out.append(s)

# ---------- degree laws ----------
def mk(d, p, label):
    d = np.asarray(d, float); p = np.asarray(p, float); p = p / p.sum(); keep = p > 1e-12
    d, p = d[keep], p[keep]; mu1 = (d * p).sum(); mu2 = (d * d * p).sum()
    return dict(d=d, p=p, mu1=mu1, mu2=mu2, nu=mu2 / mu1 - 1, label=label)
def nb(mu, cv):
    inv_r = max(cv**2 - 1 / mu, 1e-3); d = np.arange(1, int(round(mu * (1 + 8 * cv))) + 1)
    return mk(d, stats.nbinom.pmf(d, 1 / inv_r, (1 / inv_r) / (1 / inv_r + mu)), f"NB {mu}/{cv}")
def ln(mu, cv):
    s2 = np.log(1 + cv**2); m = np.log(mu) - s2 / 2; d = np.arange(1, int(round(mu * (1 + 12 * cv))) + 1)
    return mk(d, stats.lognorm.pdf(d, np.sqrt(s2), scale=np.exp(m)), f"LN {mu}/{cv}")
def pl(mu):
    d = np.arange(1, 30 * mu + 1); dmin = max(2, int(round(mu / 3)))
    f = lambda a: ((d * np.where(d >= dmin, d**-a, 0)).sum() / np.where(d >= dmin, d**-a, 0).sum()) - mu
    a = optimize.brentq(f, 1.05, 6); return mk(d, np.where(d >= dmin, d**-a, 0), f"PL {mu}")
laws = [nb(m, c) for m in (6, 12, 25, 50) for c in (0.4, 1.0, 1.4)] + [ln(m, c) for m in (6, 12, 25, 50) for c in (0.4, 1.0, 1.4)] + [pl(m) for m in (6, 12, 25)]
hostile = []
for i in range(40):                       # random finitely supported laws with nu > 1
    K = rng.integers(2, 40); d = np.unique(rng.integers(1, rng.choice([5, 30, 200]), size=K)); p = rng.dirichlet(np.ones(len(d)) * rng.choice([0.3, 1.0]))
    L = mk(d, p, f"rand{i}")
    if L["nu"] > 1.05: hostile.append(L)
hostile.append(mk([1, 50], [0.97, 0.03], "mostly leaves + hubs")); hostile.append(mk([1, 2, 3], [0.2, 0.5, 0.3], "tiny degrees")); hostile.append(mk([1, 3], [0.5, 0.5], "nu near 1.5"))
hostile = [L for L in hostile if L["nu"] > 1.0]
log(f"{len(laws)} paper laws, {len(hostile)} hostile laws (nu from {min(L['nu'] for L in hostile):.2f} to {max(L['nu'] for L in hostile):.1f})")
ALL = laws + hostile
def H(L, E): E = np.asarray(E, float); return ((L["p"] * (L["d"] - 1) / L["mu1"])[:, None] * (-np.expm1(-L["d"][:, None] * E[None, :].reshape(1, -1)))).sum(0).reshape(E.shape)
def Hp(L, E): E = np.asarray(E, float); return ((L["p"] * L["d"] * (L["d"] - 1) / L["mu1"])[:, None] * np.exp(-L["d"][:, None] * E.reshape(1, -1))).sum(0).reshape(E.shape)
def fixed_point(L): return optimize.brentq(lambda E: H(L, np.array([E]))[0] - E, 1e-9, 1.0 - 1e-12) if L["nu"] > 1 else None
def omega(L, u):
    m = int(np.ceil(14 * np.log(10) / np.log(L["nu"]))) + 5          # nu^m >= 1e14, so the linearisation error u^2/nu^m is negligible
    E = np.asarray(u, float) / L["nu"]**m
    for _ in range(m): E = H(L, E)
    return E

# ---------- T1 properties of H ----------
log("\nT1  Properties of H (claimed: increasing, strictly concave, H(0)=0, H'(0)=nu, H<1-1/mu1, H(E)<=nu E, unique fixed point E* in (0,1) with H'(E*)<1)")
worst = dict(conc=0, bound=0, Hup=0, fp=0, hp=0, nu=0); bad = []
for L in ALL:
    E = np.linspace(0, 1, 4001); h = H(L, E); hp = Hp(L, E)
    Em = np.linspace(0, min(1.0, 600.0 / L["d"].max()), 2001); hpm = Hp(L, Em)       # exp(-dE) would underflow for dE > 700, so test H' > 0 on dE <= 600
    ok = (np.all(np.diff(hp) <= 1e-12) and np.all(hpm > 0) and np.all(h <= L["nu"] * E + 1e-12) and np.all(h < 1 - 1 / L["mu1"] + 1e-12) and abs(Hp(L, np.array([0.0]))[0] - L["nu"]) < 1e-9 * L["nu"])
    Es = fixed_point(L); g = H(L, E[1:]) / E[1:]; unique = np.sum(np.diff(np.sign(g - 1)) != 0) == 1
    ok = ok and unique and Hp(L, np.array([Es]))[0] < 1 and 0 < Es < 1
    if not ok: bad.append(L["label"])
log(f"    {len(ALL) - len(bad)}/{len(ALL)} laws pass all property checks" + (f"; FAIL: {bad}" if bad else ""))

# ---------- T2 Poincare function ----------
log("\nT2  omega = phi^{-1}: H(omega(u)) = omega(nu u) globally; omega increasing; omega(u) = u + O(u^2); omega(u) -> E* as u -> inf")
res = []
for L in ALL:
    u = np.geomspace(1e-6, 1e3, 200); w = omega(L, u); Es = fixed_point(L)
    rel = np.max(np.abs(H(L, w) - omega(L, L["nu"] * u)) / np.maximum(omega(L, L["nu"] * u), 1e-300)) if True else 0
    mono = np.all(np.diff(w) >= 0); lim = abs(omega(L, np.array([1e9]))[0] - Es)
    small = omega(L, np.array([1e-4]))[0] / 1e-4 - 1          # should be O(u)
    res.append((rel, mono, lim, small, L["label"]))
log(f"    worst relative functional-equation error {max(r[0] for r in res):.1e}; monotone in all: {all(r[1] for r in res)}; worst |omega(1e6)-E*| {max(r[2] for r in res):.1e}; worst omega(1e-4)/1e-4 - 1 = {max(abs(r[3]) for r in res):.1e}")
# ---------- T3 Lemma A ----------
log("\nT3  Lemma A: sup_k |E_k - omega(nu^k delta/(nu-1))| over nu^k delta/(nu-1) in [a,b], as delta -> 0 (claimed -> 0); and the linear-phase bounds")
def orbit(L, delta, K): 
    E = 0.0; out = []
    for _ in range(K): E = delta + H(L, np.array([E]))[0]; out.append(E)
    return np.array(out)
tab = {}
for L in ALL:
    for delta in (1e-3, 1e-5, 1e-7, 1e-9):
        K = int(np.ceil(np.log((50 * (L["nu"] - 1)) / delta) / np.log(L["nu"]))) + 3; E = orbit(L, delta, K); k = np.arange(1, K + 1); x = L["nu"]**k * delta / (L["nu"] - 1)
        sel = (x >= 0.1) & (x <= 20); err = np.max(np.abs(E[sel] - omega(L, x[sel]))) if sel.any() else np.nan
        upper_ok = np.all(E <= x + 1e-15); tab.setdefault(delta, []).append((err, upper_ok, L["label"]))
for delta, v in tab.items():
    errs = [e for e, _, _ in v if not np.isnan(e)]; log(f"    delta={delta:.0e}: worst sup-error over [a,b]=[0.1,20] = {max(errs):.2e} (law {max(v, key=lambda t: (-1 if np.isnan(t[0]) else t[0]))[2]}); upper bound E_k <= delta nu^k/(nu-1) holds in all: {all(u for _, u, _ in v)}")
# rate: error should scale like delta
e3 = max(e for e, _, _ in tab[1e-3] if not np.isnan(e)); e9 = max(e for e, _, _ in tab[1e-9] if not np.isnan(e)); log(f"    error ratio delta=1e-3 vs 1e-9: {e3 / e9:.1e} (convergence to 0 confirmed if >> 1)")
# ---------- T4 Q-calligraphic and translation ----------
def Qcal(L, lam):
    kap = L["mu1"] * (L["nu"] - 1); num = den = 0.0
    for s, ps in zip(L["d"], L["p"]):
        a = omega(L, np.array([lam * s / (L["nu"] * kap)]))[0]; b = omega(L, np.array([lam * s / kap]))[0]; pi = np.exp(-L["d"] * a) - np.exp(-L["d"] * b)
        num += ps * (L["p"] * L["d"] * pi).sum(); den += ps * (L["p"] * pi).sum()
    return num / (L["mu1"] * den)
def Qtilde(L, n, K):
    num = np.zeros(K); den = np.zeros(K)
    for s, ps in zip(L["d"], L["p"]):
        E = orbit(L, s / (n * L["mu1"]), K); Em = np.concatenate([[0.0], E[:-1]]); pi = np.exp(-np.outer(L["d"], Em)) - np.exp(-np.outer(L["d"], E))
        num += ps * (L["p"] * L["d"]) @ pi; den += ps * L["p"] @ pi
    return num / (L["mu1"] * den)
log("\nT4  Q-tilde_k(n) vs Qcal(nu^k/n) on shells with nu^k/n in [1e-3, 1e3]; Qcal(lambda) as lambda -> 0 vs mu2/mu1^2; translation error")
sub = [L for L in laws if L["label"] in ("NB 12/1.0", "NB 6/0.4", "LN 25/1.4", "PL 12", "NB 50/1.4")] + hostile[:6]
for n in (1e3, 1e5, 1e7):
    worst = 0
    for L in sub:
        K = int(np.log(n * 1e3) / np.log(L["nu"])) + 2; Qt = Qtilde(L, n, K); k = np.arange(1, K + 1); lam = L["nu"]**k / n; sel = (lam > 1e-3) & (lam < 1e3)
        for kk in np.where(sel)[0]: worst = max(worst, abs(Qt[kk] - Qcal(L, lam[kk])))
    log(f"    n={n:.0e}: worst |Q-tilde - Qcal| = {worst:.1e}")
for L in sub[:6]: log(f"    {L['label']}: Qcal(1e-8) = {Qcal(L, 1e-8):.6f}  vs  mu2/mu1^2 = {L['mu2'] / L['mu1']**2:.6f}")
# ---------- T5 reversal distance: discrete and continuous, dense in n so that u(n) comes close to integers ----------
def lam_star(L):
    f = lambda t: Qcal(L, np.exp(t)) - 1; ts = np.linspace(-6, 8, 57); vals = np.array([f(t) for t in ts]); j = np.where(np.diff(np.sign(vals)) != 0)[0]
    if len(j) != 1: return None
    return np.exp(optimize.brentq(f, ts[j[0]], ts[j[0] + 1]))
def krev_cont(Qt):
    k = np.where(Qt[1:] < 1)[0]
    if len(k) == 0: return None, None
    kr = k[0] + 2; return kr, (kr - 1) + (Qt[kr - 2] - 1) / (Qt[kr - 2] - Qt[kr - 1])
log("\nT5  Discrete and continuous reversal distance vs the proposition, n dense in log scale (so that u(n) approaches integers)")
for L in [l for l in laws if l["label"] in ("NB 12/1.0", "LN 25/1.4", "PL 12", "NB 6/0.4")] + hostile[:3]:
    ls = lam_star(L)
    if ls is None: log(f"    {L['label']}: single crossing of Qcal NOT found on grid"); continue
    g = lambda t: Qcal(L, ls * L["nu"]**t); 
    def pi_fun(u):
        f = u - np.floor(u); return -f + (g(-f) - 1) / (g(-f) - g(1 - f))
    errs = []; dis = []; near = 0
    for n in np.geomspace(2e3, 2e3 * L["nu"]**3, 60):
        u = np.log(n * ls) / np.log(L["nu"]); K = int(np.ceil(u)) + 4; Qt = Qtilde(L, n, K); kr, kc = krev_cont(Qt)
        if kr is None: continue
        errs.append(abs(kc - (u + pi_fun(u)))); dist = abs(u - round(u)); dis.append((dist, kr, np.ceil(u), round(u)))
    # discrete: equals ceil(u) when dist >= 0.05, else in {m, m+1}
    far = [t for t in dis if t[0] >= 0.05]; nearl = [t for t in dis if t[0] < 0.05]
    okfar = all(t[1] == t[2] for t in far); oknear = all(t[1] in (t[3], t[3] + 1) for t in nearl)
    log(f"    {L['label']:<22} lam*={ls:8.3f} max |k~rev - (u+pi(u))| = {max(errs):.2e} (n from 2e3 to {2e3 * L['nu']**3:.1e}); discrete = ceil(u) when far from integers: {okfar} ({len(far)} n); in {{m,m+1}} when near: {oknear} ({len(nearl)} n)")
log("\nT6  Existence of lambda* (scope of the single-crossing hypothesis): number of crossings of Qcal=1 for lambda in [e^-8, e^21], all 70 laws")
cnt = {}
for L in ALL:
    lam = np.exp(np.linspace(-8, 21, 59)); Q = np.array([Qcal(L, l) for l in lam]); ok = np.isfinite(Q); Q = Q[ok]; c = int(np.sum(np.diff(np.sign(Q - 1)) != 0)); cnt[L["label"]] = (c, L["nu"], fixed_point(L))
paper = [v for k, v in cnt.items() if k in {l["label"] for l in laws}]
log(f"    27 paper laws: crossings = {sorted(set(v[0] for v in paper))}; nu range {min(v[1] for v in paper):.2f} to {max(v[1] for v in paper):.1f}; E* >= {min(v[2] for v in paper):.3f}")
none = [(k, v) for k, v in cnt.items() if v[0] == 0]; multi = [k for k, v in cnt.items() if v[0] > 1]
log(f"    hostile laws: {len(none)} with no crossing (Q stays above 1: nu in [{min(v[1] for _, v in none):.2f}, {max(v[1] for _, v in none):.2f}], E* <= {max(v[2] for _, v in none):.3f}); laws with more than one crossing: {len(multi)}")
open("../results/proof_audit.txt", "w").write("\n".join(out))
