#!/usr/bin/env python3
"""Compare the GEM strip-timing-cut A/B replays (same segments, lgem_usestriptimingcut 1 vs 0).

  python3 ana_ab.py <skim_cut.root> <skim_nocut.root> <png dir>

Both skims come from skim_trk.C. Signal = tof-L/c [30,50] minus in-time sidebands
[10,25] U [55,62]; random-cluster rates f from the out-of-time hits; eps = (m-f)/(1-f).
"""
import sys

import numpy as np
import uproot

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

KEYS = ["run", "evt", "tofc", "pred", "res", "tn", "chi", "ncl"]
SB = [(10.0, 25.0), (55.0, 62.0)]
K = 20.0 / sum(b - a for a, b in SB)
SL = ["x0", "y0", "x1", "y1"]
VARS = ["x0", "x1", "xb", "y0", "y1", "yb", "c0", "c1", "cb"]


def load(fn):
    D = uproot.open(fn)["H"].arrays(KEYS, library="np")
    t, p = D["tofc"], D["pred"]
    D["pk"] = (t >= 30) & (t < 50)
    sb = np.zeros(len(t), bool)
    for a, b in SB:
        sb |= (t >= a) & (t < b)
    D["sb"] = sb
    D["oot"] = (t < -100) | (t > 125)
    D["fid"] = (p[:, 3] < 19.9) & (p[:, 3] > -25.7) & (np.abs(p[:, 2]) < 56.4) & \
               (p[:, 1] < 19.9) & (p[:, 1] > -25.7) & (np.abs(p[:, 0]) < 56.4)
    D["ares"] = np.abs(D["res"])
    return D


def eps(D, M, sel):
    pk, sb, oot = D["pk"], D["sb"], D["oot"]
    a, b = (sel & pk).sum(), (sel & sb).sum()
    S = a - K * b
    o = ((sel & pk & M).sum() - K * (sel & sb & M).sum()) / S
    f = (sel & oot & M).sum() / max((sel & oot).sum(), 1)
    err = np.sqrt(max(o * (1 - o), 0.01) * (a + K * K * b)) / S / (1 - f)
    return (o - f) / (1 - f), err, f, S


A, B = load(sys.argv[1]), load(sys.argv[2])
png = sys.argv[3]
for lab, D in (("cut", A), ("nocut", B)):
    print(f"{lab}: hits {len(D['tofc'])}, signal {eps(D, np.ones(len(D['tofc']), bool), np.ones(len(D['tofc']), bool))[3]:.0f}, "
          f"mean clusters/slot (OOT) " + " ".join(f"{D['ncl'][D['oot'], s].mean():.1f}" for s in range(4)))

print("\n== per-slot match efficiency vs residual window (all / fiducial), cut -> nocut  [f = random rate]")
for W in (3, 5, 10):
    for s in range(4):
        row = []
        for D in (A, B):
            ALL = np.ones(len(D["tofc"]), bool)
            e1 = eps(D, D["ares"][:, s] < W, ALL)
            e2 = eps(D, D["ares"][:, s] < W, D["fid"])
            row.append(f"{e1[0]:.3f}+-{e1[1]:.3f} (f {e1[2]:.2f}) | fid {e2[0]:.3f}+-{e2[1]:.3f} (f {e2[2]:.2f})")
        print(f"  W={W:2d} {SL[s]}:  cut {row[0]}   ->   nocut {row[1]}")

print("\n== GEM1 any (x1<10 or y1<5), fiducial")
for lab, D in (("cut", A), ("nocut", B)):
    M = (D["ares"][:, 2] < 10) | (D["ares"][:, 3] < 5)
    e = eps(D, M, D["fid"])
    print(f"  {lab}: {e[0]:.3f}+-{e[1]:.3f}  f={e[2]:.2f}  S={e[3]:.0f}")

print("\n== chi2 < 50 tracked variants (eq. 8), all hits")
for lab, D in (("cut", A), ("nocut", B)):
    ALL = np.ones(len(D["tofc"]), bool)
    out = []
    for i, v in enumerate(VARS):
        M = (D["chi"][:, i] >= 0) & (D["chi"][:, i] < 50)
        e = eps(D, M, ALL)
        out.append(f"{v} {e[0]:.3f}(f{e[2]:.2f})")
    print(f"  {lab}: " + "  ".join(out))

print("\n== pair (x0 W5, x1 W10) in fiducial: eps0, eps1, eps01, none")
for lab, D in (("cut", A), ("nocut", B)):
    F = D["fid"]
    m0, m1 = D["ares"][:, 0] < 5, D["ares"][:, 2] < 10
    e0, e1 = eps(D, m0, F)[0], eps(D, m1, F)[0]
    f0, f1 = eps(D, m0, F)[2], eps(D, m1, F)[2]
    f01 = (F & D["oot"] & m0 & m1).sum() / (F & D["oot"]).sum()
    pk, sb = D["pk"], D["sb"]
    S = (F & pk).sum() - K * (F & sb).sum()
    o = ((F & pk & m0 & m1).sum() - K * (F & sb & m0 & m1).sum()) / S
    e01 = (o - e0 * f1 - e1 * f0 - (1 - e0 - e1) * f01) / (1 - f0 - f1 + f01)
    print(f"  {lab}: {e0:.3f} {e1:.3f} {e01:.3f}  none={1 - e0 - e1 + e01:.3f}")

# matched-cluster time for signal, both configs
tb = np.arange(0, 160, 4)
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
for s in range(4):
    a = ax.flat[s]
    for lab, D, col in (("cut", A, "b"), ("nocut", B, "r")):
        m = D["ares"][:, s] < 3
        hp, _ = np.histogram(D["tn"][m & D["pk"], s], tb)
        hs, _ = np.histogram(D["tn"][m & D["sb"], s], tb)
        a.step(0.5 * (tb[1:] + tb[:-1]), hp - K * hs, where="mid", color=col, label=f"{lab}: signal, |res|<3")
    for x in (55.5, 118.5):
        a.axvline(x, ls="--", color="gray")
    a.set_title(f"{SL[s]}: time of matched cluster (peak - sidebands)")
    a.set_xlabel("clust.time (ns)")
    a.legend()
fig.tight_layout()
fig.savefig(f"{png}/ab_cluster_time.png", dpi=110)
print("\n== signal matched-cluster time fractions (|res|<3), nocut")
for s in range(4):
    D = B
    m = D["ares"][:, s] < 3
    tt = D["tn"][:, s]
    def S(x):
        return (x & D["pk"]).sum() - K * (x & D["sb"]).sum()
    tot = S(m)
    print(f"  {SL[s]}: t<55.5 {S(m & (tt < 55.5)) / tot:.3f}  55.5-118.5 {S(m & (tt >= 55.5) & (tt <= 118.5)) / tot:.3f}  >118.5 {S(m & (tt > 118.5)) / tot:.3f}")
