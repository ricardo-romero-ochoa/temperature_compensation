#!/usr/bin/env python3
"""Figure 4: finite-equilibration law CV^2=(1/m)(1+2 eps K/(1+K)) and the exact
activity identity (CV^2_step-1) K = 4K^2/(1+K)^2. Reproduces finite.png."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tempcomp.oscillators import gated_chain_cv2_formula, activity_identity_rhs

OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "fig4_finite.png")


def main():
    K = 20.0
    eps = np.geomspace(1e-3, 1, 200)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for m, c in [(1, "#c0392b"), (3, "#e67e22"), (6, "#16a085")]:
        ax[0].plot(eps, gated_chain_cv2_formula(m, eps, K), color=c, lw=2, label=f"m={m}")
        ax[0].axhline(1 / m, ls=":", color=c, lw=1)
    ax[0].set_xscale("log")
    ax[0].set_xlabel(r"equilibration parameter $\epsilon=p/(f+b)$")
    ax[0].set_ylabel(r"CV$^2_{\mathrm{chain}}$"); ax[0].legend()
    ax[0].text(0.02, 0.97, "A", transform=ax[0].transAxes, fontsize=15,
               fontweight="bold", va="top")

    for K, c in [(2.0, "#2c3e50"), (20.0, "#8e44ad")]:
        act = np.geomspace(2, 300, 100)
        ax[1].plot(act, activity_identity_rhs(K) / act, color=c, lw=2,
                   label=rf"$K_{{\mathrm{{eq}}}}={K:.0f}$")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel(r"dynamical activity per step $\mathcal{A}_{\mathrm{step}}$")
    ax[1].set_ylabel(r"excess variance  CV$^2_{\mathrm{step}}-1$")
    ax[1].legend()
    ax[1].text(0.02, 0.97, "B", transform=ax[1].transAxes, fontsize=15,
               fontweight="bold", va="top")
    plt.tight_layout(); plt.savefig(OUT, dpi=300)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
