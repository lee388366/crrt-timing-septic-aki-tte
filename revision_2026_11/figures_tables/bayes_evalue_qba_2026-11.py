# Bayesian reanalysis (7 priors), E-values and quantitative bias analysis from summary estimates only (unrounded HR and BCa 95% CI)
from math import log, exp, sqrt
from statistics import NormalDist
N = NormalDist()
PRIORS = [("Data-dominant","Primary",0,1),("Neutral","Primary",0,0.355),("Sceptical","Primary",0,0.205),
          ("Mildly optimistic","Sensitivity",-0.105,0.205),("Mildly pessimistic","Sensitivity",0.095,0.205),
          ("Discounted-historical neutral","Sensitivity",0,0.25),("ELAIN-informed stress test","Sensitivity",-0.223,0.3)]
EST = {"primary": (1.03083625505799, 0.885126347369368, 1.20246385510249), "full": (1.01247832042097, 0.900632494698885, 1.15263105169048)}   # primary analysis (weights truncated at the 99th percentile)
rr = lambda h: (1 - 0.5 ** sqrt(h)) / (1 - 0.5 ** sqrt(1 / h))            # common-outcome HR→RR
ev = lambda r: (1 / r if r < 1 else r) + sqrt((1 / r if r < 1 else r) * ((1 / r if r < 1 else r) - 1))
for pop, (hr, lo, hi) in EST.items():
    b = log(hr); se = (log(hi) - log(lo)) / 3.92
    print(f"[{pop}] logHR {b:.4f}  SE {se:.4f}  E-value point {ev(rr(hr)):.2f}  CI {1.0 if lo < 1 < hi else ev(rr(hi if hr < 1 else lo)):.2f}")
    for nm, role, m, s in PRIORS:
        v = 1 / (1 / s**2 + 1 / se**2); mu = v * (m / s**2 + b / se**2); sd = sqrt(v)
        ben = N.cdf((log(.9) - mu) / sd); harm = 1 - N.cdf((log(1.1) - mu) / sd)
        print(f"   {nm:30s} {role:11s} {m:+.3f} ({s:.3f}) | median HR {exp(mu):.3f} | benefit {100*ben:.1f} | ROPE {100*(1-ben-harm):.1f} | harm {100*harm:.1f}")
    if pop == "primary":
        r = rr(hr)
        for t in (0.9, 1.1):
            B = max(r / rr(t), rr(t) / r); print(f"   QBA equal-strength parameter to reach HR {t}: {B + sqrt(B*B - B):.2f}")
