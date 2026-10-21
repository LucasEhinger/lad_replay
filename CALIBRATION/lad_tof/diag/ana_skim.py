#!/usr/bin/env python3
"""First look at the skim_trk.C output: sideband-subtracted proton-peak quantities.

  python3 ana_skim.py <skim dir> <png dir>

Peak  = tof - L/c in [30, 50] ns.
Bkg   = in-time sidebands [0, 25] U [55, 80] ns (same GEM-time acceptance as the
        peak, so accidental real tracks and noise are both subtracted), scaled by
        20/50.
All "signal" numbers below are peak - k * sideband.
"""
import glob
import sys

import numpy as np
import uproot

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SLOT = ["x0", "y0", "x1", "y1"]
VARS = ["x0", "x1", "xb", "y0", "y1", "yb", "c0", "c1", "cb"]
PK = (30.0, 50.0)
SB = [(0.0, 25.0), (55.0, 80.0)]
K = (PK[1] - PK[0]) / sum(b - a for a, b in SB)

skdir, pngdir = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(f"{skdir}/skim_*.root"))
arrs = [uproot.open(f)["H"].arrays(library="np") for f in files]
D = {k: np.concatenate([a[k] for a in arrs]) for k in arrs[0]}
n = len(D["tofc"])
t = D["tofc"]
pk = (t >= PK[0]) & (t < PK[1])
sb = np.zeros(n, bool)
for a, b in SB:
    sb |= (t >= a) & (t < b)
print(f"hits {n}, peak {pk.sum()}, sidebands {sb.sum()} (k = {K:.3f}), signal {pk.sum() - K * sb.sum():.0f}")


def S(mask):
    """Background-subtracted count and its error."""
    a, b = (mask & pk).sum(), (mask & sb).sum()
    return a - K * b, np.sqrt(a + K * K * b)


def frac(num, den=None):
    den = np.ones(n, bool) if den is None else den
    s1, _ = S(num & den)
    s0, _ = S(den)
    return s1 / s0 if s0 > 0 else np.nan


def shist(x, bins, mask=None):
    m = np.ones(n, bool) if mask is None else mask
    hp, _ = np.histogram(x[m & pk], bins)
    hs, _ = np.histogram(x[m & sb], bins)
    return hp - K * hs, np.sqrt(hp + K * K * hs), hp, K * hs


res = D["res"]
chi = D["chi"]
trk = (chi >= 0) & (chi < 50)

print("\n== tracking fraction (sideband-subtracted, no noise correction; chi2 < 50)")
for i, v in enumerate(VARS):
    print(f"  {v}: {frac(trk[:, i]):.3f}")

print("\n== nearest-cluster |residual| (no road) for signal hits")
print("  slot   none   <1cm   <2cm   <3cm   <5cm  <10cm  <30cm   n_clusters(mean, peak)")
for s in range(4):
    r = np.abs(res[:, s])
    row = [frac(r > 900)] + [frac(r < c) for c in (1, 2, 3, 5, 10, 30)]
    print(f"  {SLOT[s]:4s} " + " ".join(f"{x:6.3f}" for x in row) + f"   {D['ncl'][pk, s].mean():.1f}")

# random-cluster expectation: fraction of OOT hits with a cluster within c
oot = (t < -100) | (t > 125)
print("\n== same for out-of-time hits (pure random-cluster rate)")
for s in range(4):
    r = np.abs(res[oot, s])
    print(f"  {SLOT[s]:4s} " + " ".join(f"{(r < c).mean():6.3f}" for c in (1, 2, 3, 5, 10, 30)))

# residual distributions
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
bins = np.linspace(-30, 30, 121)
for s in range(4):
    h, e, hp, hb = shist(res[:, s], bins)
    a = ax.flat[s]
    c = 0.5 * (bins[1:] + bins[:-1])
    a.step(c, hp, where="mid", label="peak", color="k")
    a.step(c, hb, where="mid", label="sideband x k", color="gray")
    a.errorbar(c, h, e, fmt=".", color="r", label="peak - bkg")
    a.set_title(f"{SLOT[s]}: nearest cluster - prediction")
    a.set_xlabel("cm")
    a.legend()
fig.tight_layout()
fig.savefig(f"{pngdir}/resid_nearest.png", dpi=110)

# cluster time of matched clusters
print("\n== time of the nearest cluster (|res| < 3 cm), signal hits")
tb = np.arange(36, 136, 4)
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
for s in range(4):
    m = np.abs(res[:, s]) < 3
    h, e, hp, hb = shist(D["tn"][:, s], tb, m)
    a = ax.flat[s]
    c = 0.5 * (tb[1:] + tb[:-1])
    a.step(c, hp, where="mid", color="k", label="peak")
    a.step(c, hb, where="mid", color="gray", label="sideband x k")
    a.errorbar(c, h, e, fmt=".", color="r", label="peak - bkg")
    for x in (55.5, 118.5):
        a.axvline(x, ls="--", color="b")
    a.set_title(f"{SLOT[s]}: cluster time, |res| < 3 cm")
    a.set_xlabel("clust.time (ns)")
    a.legend()
    print(f"  {SLOT[s]}: " + " ".join(f"{int(lo)}:{v:.0f}" for lo, v in zip(tb[:-1], h)))
fig.tight_layout()
fig.savefig(f"{pngdir}/cluster_time_matched.png", dpi=110)

# cluster time vs tof (all hits, matched within 2 cm)
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for j, s in enumerate((0, 2)):
    m = np.abs(res[:, s]) < 2
    ax[j].hist2d(t[m], D["tn"][m, s], bins=[np.arange(-150, 176, 5), np.arange(36, 136, 2)], cmap="viridis")
    ax[j].set_xlabel("tof - L/c (ns)")
    ax[j].set_ylabel("clust.time (ns)")
    ax[j].set_title(f"{SLOT[s]}: nearest cluster within 2 cm")
fig.tight_layout()
fig.savefig(f"{pngdir}/cluster_time_vs_tof.png", dpi=110)

# conditional: match within 3 cm
print("\n== P(row matched within 3 cm | column matched within 3 cm), signal hits")
mt = np.abs(res) < 3
print("        " + " ".join(f"{s:>6s}" for s in SLOT))
for a_ in range(4):
    print(f"  {SLOT[a_]:4s}  " + " ".join(f"{frac(mt[:, a_], mt[:, b]):6.3f}" for b in range(4)) +
          f"   single {frac(mt[:, a_]):.3f}")

np.savez_compressed(f"{pngdir}/skim_cache.npz", **D)
