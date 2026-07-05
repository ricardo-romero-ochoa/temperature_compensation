#!/usr/bin/env python3
"""Figure 2: single-timescale closed form -- compensated coherence vs antagonist
activation gap G, and the N_coh(chi) frontier. Reproduces theorem.png."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tempcomp.oscillators import semimarkov_coherence, compensated_delta

OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "fig2_closedform.png")
N, E_A = 4, 4.0


def chi_of_delta(N, E_A, G, Delta, T0=1.0):
    tau = N + Delta
    return abs(N * E_A - G * Delta) / (T0 ** 2 * tau)


def main():
    Gs = np.linspace(0.3, 12, 200)
    Ncomp = semimarkov_coherence(N, compensated_delta(N, E_A, Gs), 0.0)
    ceiling, floor = N / np.pi, 1 / np.pi

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].plot(Gs, Ncomp, "-", color="#2c3e50", lw=2)
    ax[0].axhline(ceiling, ls="--", color="gray")
    ax[0].axhline(floor, ls=":", color="red")
    ax[0].text(8, ceiling + 0.03, r"ceiling $N/\pi$", color="gray", fontsize=9)
    ax[0].text(6, floor + 0.03, r"floor $1/\pi$", color="red", fontsize=9)
    ax[0].set_xlabel(r"antagonist activation gap $G=E_B-E_A-E_e$")
    ax[0].set_ylabel(r"compensated coherence $\mathcal{N}^*$")
    ax[0].set_ylim(0, ceiling + 0.15)
    ax[0].text(0.02, 0.97, "A", transform=ax[0].transAxes, fontsize=15,
               fontweight="bold", va="top")

    for G, lab, col in [(2.0, "G=2 weak", "#c0392b"),
                        (4.0, "G=4", "#e67e22"),
                        (8.0, "G=8 strong", "#27ae60")]:
        Ds = np.linspace(0, 14, 400)
        chis = np.array([chi_of_delta(N, E_A, G, D) for D in Ds])
        ncs = semimarkov_coherence(N, Ds, 0.0)
        ax[1].plot(chis, ncs, "-", color=col, lw=2, label=lab)
        Dc = compensated_delta(N, E_A, G)
        ax[1].plot(0, semimarkov_coherence(N, Dc, 0.0), "o", color=col, ms=7)
    ax[1].axhline(floor, ls=":", color="red")
    ax[1].set_xlim(-0.2, 6); ax[1].legend()
    ax[1].set_xlabel(r"compensation deficit $\chi$")
    ax[1].set_ylabel(r"coherence $\mathcal{N}$")
    ax[1].text(0.02, 0.97, "B", transform=ax[1].transAxes, fontsize=15,
               fontweight="bold", va="top")
    plt.tight_layout(); plt.savefig(OUT, dpi=300)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
