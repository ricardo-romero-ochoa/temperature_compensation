#!/usr/bin/env python3
"""Figure 1: minimal two-loop oscillator -- compensation vs coherence trade-off,
achievable region, and activity-vs-dissipation. Reproduces final.png."""
import os, sys, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from tempcomp.markov import generator, eigen_coherence, entropy_production, dynamical_activity, steady_state
from tempcomp.oscillators import two_loop_edges, N_STATES_TWO_LOOP
from tempcomp.analysis import compensation_deficit_omega

OUT = os.path.join(os.path.dirname(__file__), "..", "figures", "fig1_twoloop.png")
T0 = 1.0


def observe(T, g, **kw):
    e = two_loop_edges(T, g, **kw)
    L = generator(e, N_STATES_TWO_LOOP)
    pi = steady_state(L)
    omega, coh = eigen_coherence(L)
    return dict(omega=omega, coh=coh,
                sigma=entropy_production(e, pi),
                K=dynamical_activity(e, pi))


def main():
    # (1) single-knob scan
    gs = np.geomspace(0.2, 60, 60)
    rows = []
    for g in gs:
        ob = observe(T0, g)
        chi = compensation_deficit_omega(two_loop_edges, T0, g=g)
        rows.append((g, ob["omega"], ob["coh"], abs(chi), ob["sigma"], ob["K"]))
    g_, w_, coh_, chi_, sig_, K_ = np.array(rows).T
    gstar = g_[np.argmin(chi_)]

    # (2) achievable region over parameters
    pts = []
    for bias, g, ep, intr in itertools.product(
            [0.01, 0.02, 0.04, 0.08], np.geomspace(0.1, 120, 26),
            [74, 200, 600], [0.3, 1.0, 3.0]):
        ob = observe(T0, g, bias=bias, entry_pref=ep, internal=intr)
        chi = compensation_deficit_omega(two_loop_edges, T0, g=g, bias=bias,
                                         entry_pref=ep, internal=intr)
        if ob["coh"] > 1e-3 and np.isfinite(chi):
            pts.append((abs(chi), ob["coh"], ob["sigma"]))
    pts = np.array(pts)
    chi_all, coh_all, sig_all = pts.T
    bins = np.linspace(0, 4, 21)
    env_x, env_y = [], []
    for i in range(len(bins) - 1):
        m = (chi_all >= bins[i]) & (chi_all < bins[i + 1])
        if m.sum():
            env_x.append(0.5 * (bins[i] + bins[i + 1]))
            env_y.append(coh_all[m].max())

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))
    axL, axR = ax[0], ax[0].twinx()
    axL.plot(g_, np.minimum(chi_, 8), "o-", color="#c0392b", ms=3)
    axR.plot(g_, coh_, "s-", color="#2c3e50", ms=3)
    axL.axvline(gstar, ls="--", color="gray", lw=1)
    axL.set_xscale("log"); axL.set_ylim(0, 8)
    axL.set_xlabel("antagonistic-loop competitiveness  g")
    axL.set_ylabel(r"$|d\ln\omega/dT|$ (capped)", color="#c0392b")
    axR.set_ylabel(r"coherence $\mathcal{N}$", color="#2c3e50")
    axL.text(0.02, 0.97, "A", transform=axL.transAxes, fontsize=15,
             fontweight="bold", va="top")

    m = chi_all <= 5
    sc = ax[1].scatter(chi_all[m], coh_all[m], c=np.log10(sig_all[m]),
                       cmap="viridis", s=14, alpha=0.55)
    ax[1].plot(env_x, env_y, "r-o", ms=4, lw=2, label=r"max achievable $\mathcal{N}$")
    ax[1].set_xlim(0, 5); ax[1].legend(loc="lower right", fontsize=9)
    ax[1].set_xlabel(r"compensation deficit $\chi$")
    ax[1].set_ylabel(r"coherence $\mathcal{N}$")
    ax[1].text(0.02, 0.97, "B", transform=ax[1].transAxes, fontsize=15,
               fontweight="bold", va="top")
    plt.colorbar(sc, ax=ax[1]).set_label(r"$\log_{10}\sigma$")

    ax[2].plot(g_, sig_ / sig_[0], "o-", color="#8e44ad", ms=3, label=r"$\sigma/\sigma_0$")
    ax[2].plot(g_, K_ / K_[0], "^-", color="#16a085", ms=3, label=r"$\mathcal{K}/\mathcal{K}_0$")
    ax[2].axvline(gstar, ls="--", color="gray", lw=1)
    ax[2].set_xscale("log"); ax[2].legend(fontsize=9)
    ax[2].set_xlabel("antagonistic-loop competitiveness  g")
    ax[2].set_ylabel("cost relative to uncompensated")
    ax[2].text(0.02, 0.97, "C", transform=ax[2].transAxes, fontsize=15,
               fontweight="bold", va="top")
    plt.tight_layout(); plt.savefig(OUT, dpi=300)
    print("wrote", OUT)
    print(f"compensation g*={gstar:.2f}: chi {chi_[0]:.2f}->{chi_.min():.3f}, "
          f"coh {coh_[0]:.3f}->{coh_[np.argmin(chi_)]:.3f}, "
          f"sigma x{sig_[np.argmin(chi_)]/sig_[0]:.2f}")


if __name__ == "__main__":
    main()
