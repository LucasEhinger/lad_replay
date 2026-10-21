#!/usr/bin/env python3
"""Coincidence-proton ypos distributions per paddle, with and without a GEM track.

  python3 plot_paddle_dists.py <skim cache .npz> <out dir> [variant ...]

For each back plane (001, 101) and paddle: the ypos distribution of the
coincidence protons (tof-L/c [30,50] minus the in-time sidebands [10,25] U
[55,62]), split into
  with track    = hits with a real track: (tracked - f * all) / (1 - f), with the
                  GEM-noise fraction f from the out-of-time hits of that paddle,
  without track = all - with track.
Variants are chi2 < 50 tracks: x0 x1 xb y0 y1 yb c0 c1 cb (default: x0 x1 cb).
The vertical lines mark where the line from the central foil leaves the GEM0
(orange) and GEM1 (blue) active area at the paddle centre (plot_gem_acceptance.py).
"""
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/work/submit/ehingerl/hallc/analysis/written_docs/background_subtraction")
from plot_gem_acceptance import accept, lab_point  # noqa: E402

VARS = ["x0", "x1", "xb", "y0", "y1", "yb", "c0", "c1", "cb"]
NAMES = {"x0": "1D_x_GEM0", "x1": "1D_x_GEM1", "xb": "1D_x_GEMboth", "y0": "1D_y_GEM0", "y1": "1D_y_GEM1",
         "yb": "1D_y_GEMboth", "c0": "1D_GEM0", "c1": "1D_GEM1", "cb": "1D_GEMboth"}
PLANE_NAME = {1: "001", 3: "101"}
C_ALL, C_TRK, C_NOT = "#1a1a19", "#2a78d6", "#eb6834"


def edges_on_paddle(plane, paddle, layer):
    """ypos range inside GEM<layer> at the paddle centre, central foil."""
    y = np.linspace(-230, 230, 4601)
    P = lab_point(plane, np.full_like(y, 110.0 - 22.0 * paddle), y)
    a, _ = accept(layer, np.array([0.0, 0.0, 0.0]), P)
    return (y[a].min(), y[a].max()) if a.any() else None


def main():
    npz, out = sys.argv[1], sys.argv[2]
    variants = sys.argv[3:] or ["x0", "x1", "cb"]
    D = dict(np.load(npz))
    t, pl, pd, yp, chi = D["tofc"], D["pl1"], D["pd1"], D["y1"], D["chi"]
    pk = (t >= 30) & (t < 50)
    sb = ((t >= 10) & (t < 25)) | ((t >= 55) & (t < 62))
    K = 20.0 / 22.0
    oot = (t < -100) | (t > 125)
    bins = np.arange(-225, 226, 30.0)
    c = 0.5 * (bins[1:] + bins[:-1])
    for v in variants:
        trk = (chi[:, VARS.index(v)] >= 0) & (chi[:, VARS.index(v)] < 50)
        fig, axes = plt.subplots(2, 11, figsize=(22, 6.5), sharex=True)
        summary = []
        for r, plane in enumerate((1, 3)):
            ymax = 0
            for p in range(11):
                ax = axes[r, p]
                sel = (pl == plane) & (pd == p)
                f = (sel & oot & trk).sum() / max((sel & oot).sum(), 1)

                def sig(m):
                    hp, _ = np.histogram(yp[sel & pk & m], bins)
                    hs, _ = np.histogram(yp[sel & sb & m], bins)
                    return hp - K * hs, hp + K * K * hs

                a, va = sig(np.ones(len(t), bool))
                tr, vt = sig(trk)
                real = (tr - f * a) / (1 - f)
                err_real = np.sqrt(np.maximum(vt + f * f * va, 0)) / (1 - f)
                notrk = a - real
                ax.step(c, a, where="mid", color=C_ALL, lw=1.2, label="all")
                ax.errorbar(c, real, err_real, fmt="o", ms=2.5, color=C_TRK, lw=0.8, label="with track")
                ax.step(c, notrk, where="mid", color=C_NOT, lw=1.2, label="without track")
                ax.axhline(0, color="#c3c2b7", lw=0.6)
                for layer, col in ((0, C_NOT), (1, C_TRK)):
                    e = edges_on_paddle(PLANE_NAME[plane], p, layer)
                    if e:
                        for x in e:
                            ax.axvline(x, color=col, ls="--", lw=0.8)
                S, Sr = a.sum(), real.sum()
                summary.append((PLANE_NAME[plane], p, S, Sr / S if S > 0 else np.nan, f))
                ax.set_title(f"{PLANE_NAME[plane]} pd {p}\nS={S:.0f}, eff={Sr / S:.2f}" if S > 0 else
                             f"{PLANE_NAME[plane]} pd {p}", fontsize=8)
                ax.tick_params(labelsize=7)
                ymax = max(ymax, a.max())
            for p in range(11):
                axes[r, p].set_ylim(-0.25 * ymax, 1.15 * ymax)
        for ax in axes[1]:
            ax.set_xlabel("ypos (cm)", fontsize=8)
        axes[0, 0].set_ylabel("plane 001: signal / 30 cm", fontsize=8)
        axes[1, 0].set_ylabel("plane 101: signal / 30 cm", fontsize=8)
        h, l = axes[0, 0].get_legend_handles_labels()
        h += [plt.Line2D([], [], color=C_NOT, ls="--"), plt.Line2D([], [], color=C_TRK, ls="--")]
        l += ["GEM0 edge (z=0)", "GEM1 edge (z=0)"]
        fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=5, fontsize=9, frameon=False)
        fig.suptitle(f"Coincidence protons per paddle, {NAMES[v]} chi2 < 50 (noise-corrected), 10sigma_new", y=0.995,
                     fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.9))
        fig.savefig(f"{out}/paddle_ypos_{v}.png", dpi=100)
        fig.savefig(f"{out}/paddle_ypos_{v}.pdf")
        print(f"\n{NAMES[v]}: plane paddle signal eff f")
        for row in summary:
            print("  %s %2d %6.0f %5.2f %4.2f" % row)


if __name__ == "__main__":
    main()
