#!/usr/bin/env python3
"""Figure 3: no precision penalty -- CV^2_min = 1/m is flat in compensation
strength R up to the activation-budget ceiling. Reproduces nofloor.png.

Also verifies numerically (via the exact gated-chain sub-generator) that CV^2
stays ~1/m as R is swept by the trap enthalpy dH."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tempcomp import fpt
from tempcomp.oscillators import gated_chain_subgen
from tempcomp.analysis import dlnmu_dT

OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "fig3_nofloor.png")
Emax, Ep, T = 12.0, 1.0, 1.0


def main():
    # numerical verification: CV^2 vs R (swept by dH) for several m
    print("Verification: CV^2 vs required sensitivity R (gated-chain, exact FPT)")
    for m in [1, 2, 4]:
        line = []
        for dH in [2, 4, 6, 8, 10, 12]:
            subgen = lambda TT: gated_chain_subgen(TT, m, p0=1.0, E_p=0.5,
                                                   B=3000.0, K1=60.0, dH=dH)
            cv = fpt.cv2(subgen(T))
            R = dlnmu_dT(subgen, T)
            line.append((R, cv))
        rs = ", ".join(f"R={r:.2f}:CV2={c:.3f}" for r, c in line)
        print(f"  m={m} (1/m={1/m:.3f}): {rs}")

    Rmax = (Emax - Ep) / T ** 2
    plt.figure(figsize=(7.4, 5))
    for m, c in [(1, "#c0392b"), (2, "#e67e22"), (4, "#16a085"), (8, "#2c3e50")]:
        plt.hlines(1 / m, 0, Rmax, color=c, lw=2.5)
        plt.text(Rmax + 0.25, 1 / m, rf"$m={m}$", va="center", fontsize=10, color=c)
    plt.axvline(Rmax, ls="--", color="gray", lw=1)
    plt.annotate(r"$R_{\max}$", xy=(Rmax, 0.05), xytext=(Rmax - 1.6, 0.05),
                 fontsize=10, color="gray", va="center",
                 arrowprops=dict(arrowstyle="->", color="gray"))
    plt.fill_betweenx([0, 1.08], Rmax, Rmax + 2.5, color="0.85", alpha=0.6)
    plt.text(Rmax + 1.25, 0.9, "infeasible", ha="center", fontsize=9, color="0.4")
    plt.xlabel(r"required temperature sensitivity $R=d\ln\mu_0/dT$")
    plt.ylabel(r"min CV$^2$ of compensating element")
    plt.xlim(0, Rmax + 2.5); plt.ylim(0, 1.08)
    plt.tight_layout(); plt.savefig(OUT, dpi=300)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
