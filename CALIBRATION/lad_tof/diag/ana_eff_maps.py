#!/usr/bin/env python3
"""Matching efficiency per GEM slot versus hit/geometry variables, from the skim cache.

  python3 ana_eff_maps.py <png dir>      (reads <png dir>/skim_cache.npz from ana_skim.py)

A slot "matches" a hit when its nearest cluster is within W cm of the vertex -> hodo
prediction (no chi2). In any bin:
  m = sideband-subtracted matched fraction of the peak signal,
  f = matched fraction of out-of-time hits in the same bin (random clusters only),
  eps = (m - f) / (1 - f)                                   (eq. 8 with a residual window)
"""
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

png = sys.argv[1]
D = dict(np.load(f"{png}/skim_cache.npz"))
SLOT = ["x0", "y0", "x1", "y1"]
t = D["tofc"]
n = len(t)
pk = (t >= 30) & (t < 50)
sb = ((t >= 0) & (t < 25)) | ((t >= 55) & (t < 80))
K = 20.0 / 50.0
oot = (t < -100) | (t > 125)
res = np.abs(D["res"])


def eps(match, sel):
    a, b = (sel & pk).sum(), (sel & sb).sum()
    S = a - K * b
    if S <= 20:
        return np.nan, np.nan, S
    ma, mb = (sel & pk & match).sum(), (sel & sb & match).sum()
    m = (ma - K * mb) / S
    no = (sel & oot).sum()
    f = (sel & oot & match).sum() / no if no > 20 else np.nan
    e = (m - f) / (1 - f)
    # error: binomial-ish on m with the background inflation, f error ignored
    err = np.sqrt(max(m * (1 - m), 0.01) * (a + K * K * b)) / S / (1 - f)
    return e, err, S


def scan(var, edges, W=10.0, label=None, sel=None, fname=None):
    sel = np.ones(n, bool) if sel is None else sel
    c = 0.5 * (edges[1:] + edges[:-1])
    out = {}
    for s in range(4):
        r = [eps(res[:, s] < W, sel & (var >= lo) & (var < hi)) for lo, hi in zip(edges[:-1], edges[1:])]
        out[SLOT[s]] = np.array(r)
    print(f"\n== eps(W={W:g} cm) vs {label}")
    print(f"  {'bin':>12s} {'S':>7s} " + " ".join(f"{s:>12s}" for s in SLOT))
    for i in range(len(c)):
        print(f"  {edges[i]:5.0f}-{edges[i+1]:<6.0f} {out['x0'][i, 2]:7.0f} " +
              " ".join(f"{out[s][i, 0]:6.2f}+-{out[s][i, 1]:4.2f}" for s in SLOT))
    if fname:
        fig, ax = plt.subplots(figsize=(8, 5))
        for s, col in zip(SLOT, ["C0", "C1", "C2", "C3"]):
            ax.errorbar(c, out[s][:, 0], out[s][:, 1], fmt="o-", color=col, label=s, ms=4)
        ax.set_ylim(0, 1.2)
        ax.axhline(1, color="gray", lw=0.5)
        ax.set_xlabel(label)
        ax.set_ylabel(f"match efficiency, |res| < {W:g} cm")
        ax.legend()
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(f"{png}/{fname}.png", dpi=110)
    return out


print("== overall eps vs residual window W")
for W in (1, 2, 3, 5, 7, 10, 15, 20, 30):
    print(f"  W={W:3d}: " + "  ".join(f"{SLOT[s]} {eps(res[:, s] < W, np.ones(n, bool))[0]:.3f}" for s in range(4)))

scan(D["pred"][:, 2], np.arange(-60, 61, 10.0), label="GEM1 predicted V (x) position (cm)", fname="eff_vs_predx1")
scan(D["pred"][:, 3], np.arange(-35, 31, 5.0), label="GEM1 predicted U (y) position (cm)", fname="eff_vs_predy1")
scan(D["pred"][:, 0], np.arange(-60, 61, 10.0), label="GEM0 predicted V (x) position (cm)", fname="eff_vs_predx0")
scan(D["pred"][:, 1], np.arange(-35, 31, 5.0), label="GEM0 predicted U (y) position (cm)", fname="eff_vs_predy0")
scan(D["y1"], np.arange(-225, 226, 25.0), label="hodo ypos (back plane, cm)", fname="eff_vs_ypos")
scan(D["pd1"].astype(float), np.arange(-0.5, 11, 1.0), label="paddle (back plane)", fname="eff_vs_paddle")
scan(D["pl1"].astype(float), np.array([0.5, 1.5, 2.5, 3.5]), label="back plane (1=001, 3=101)")
scan(t, np.arange(30, 51, 4.0), label="tof - L/c (ns)", fname="eff_vs_tof")
scan(D["vz"], np.array([-20, -15, -5, 5, 15, 20.0]), label="vertex z (cm)")
runs = np.unique(D["run"])
redges = np.concatenate([runs[::8].astype(float) - 0.5, [runs[-1] + 0.5]])
scan(D["run"].astype(float), redges, label="run", fname="eff_vs_run")
scan(D["ncl"][:, 2].astype(float), np.array([0, 4, 7, 10, 14, 20, 100.0]), label="GEM1 V cluster count")
scan(D["ed1"], np.array([0, 10, 20, 30, 45, 60, 100, 300.0]), label="back-plane edep (MeV)")
scan(D["ed0"], np.array([0, 100, 150, 200, 300, 1000.0]), label="front-plane edep (MeV)")
scan(D["nprot"].astype(float), np.array([0.5, 1.5, 2.5, 4.5, 100]), label="proton hits in event")

# 2D maps on GEM1 and GEM0 (W = 10)
for L, (sx, sy) in enumerate([(0, 1), (2, 3)]):
    xe, ye = np.arange(-60, 61, 15.0), np.arange(-35, 31, 7.5)
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    for j, s in enumerate((sx, sy)):
        M = np.full((len(xe) - 1, len(ye) - 1), np.nan)
        for i in range(len(xe) - 1):
            for k in range(len(ye) - 1):
                sel = (D["pred"][:, sx] >= xe[i]) & (D["pred"][:, sx] < xe[i + 1]) & \
                      (D["pred"][:, sy] >= ye[k]) & (D["pred"][:, sy] < ye[k + 1])
                e, err, S = eps(res[:, s] < 10, sel)
                if S > 60 and err < 0.25:
                    M[i, k] = e
        im = ax[j].pcolormesh(xe, ye, M.T, vmin=0, vmax=1.1, cmap="viridis")
        ax[j].set_title(f"GEM{L} {SLOT[s]} match eff (|res|<10)")
        ax[j].set_xlabel("pred V (x) [cm]")
        ax[j].set_ylabel("pred U (y) [cm]")
        fig.colorbar(im, ax=ax[j])
    a, b = pk, sb
    H1, _, _ = np.histogram2d(D["pred"][a, sx], D["pred"][a, sy], [xe, ye])
    H2, _, _ = np.histogram2d(D["pred"][b, sx], D["pred"][b, sy], [xe, ye])
    im = ax[2].pcolormesh(xe, ye, (H1 - K * H2).T, cmap="magma")
    ax[2].set_title(f"GEM{L}: signal proton hits (peak - sb)")
    fig.colorbar(im, ax=ax[2])
    fig.tight_layout()
    fig.savefig(f"{png}/effmap_gem{L}.png", dpi=110)
