"""
oscillators.py -- the specific models and constructions of the paper.

Arrhenius rates are k = A * exp(-E / T), with k_B = 1.
"""
from __future__ import annotations
import numpy as np


def arrhenius(A, E, T):
    return A * np.exp(-E / T)


# ---------------------------------------------------------------------------
# 1) Minimal two-loop oscillator  (Sec. 4, Fig. 1)
#    Productive 4-ring 0->1->2->3->0 + entropic side-excursion at node 0:
#    entry 0->4 (activation E_B), fast internal return 4->5->0 (activation E_A).
# ---------------------------------------------------------------------------
def two_loop_edges(T, g, bias=0.03, entry_pref=200.0, internal=0.3,
                   E_A=4.0, E_B=9.0):
    """Edge list [(src, dst, rate), ...] for the 6-state two-loop model."""
    e = []
    ring = [(0, 1), (1, 2), (2, 3), (3, 0)]
    for s, d in ring:
        e.append((s, d, arrhenius(1.0, E_A, T)))
        e.append((d, s, arrhenius(bias, E_A, T)))
    # side excursion sharing node 0
    e.append((0, 4, g * arrhenius(entry_pref, E_B, T)))
    e.append((4, 0, g * arrhenius(entry_pref * bias, E_B, T)))
    e.append((4, 5, arrhenius(internal, E_A, T)))
    e.append((5, 4, arrhenius(internal * bias, E_A, T)))
    e.append((5, 0, arrhenius(internal, E_A, T)))
    e.append((0, 5, arrhenius(internal * bias, E_A, T)))
    return e


N_STATES_TWO_LOOP = 6


# ---------------------------------------------------------------------------
# 2) Semi-Markov (single-trap) reduction closed forms  (Sec. 4, Eq. 6-9)
#    Units a = k = 1. Delta = mean delay, mu_e = excursion mean time.
# ---------------------------------------------------------------------------
def semimarkov_tau_var(N, Delta, mu_e):
    """Closed-form cycle mean and variance for the single-trap reduction."""
    tau = N + Delta
    var = N + 2.0 * Delta + 2.0 * mu_e * Delta + Delta ** 2
    return tau, var


def semimarkov_coherence(N, Delta, mu_e):
    """Renewal coherence N_coh = tau^2 / (pi Var)."""
    tau, var = semimarkov_tau_var(N, Delta, mu_e)
    return tau ** 2 / (np.pi * var)


def compensated_delta(N, E_A, G):
    """Delta* = N E_A / G, G = E_B - E_A - E_e (compensation condition)."""
    return N * E_A / G


# ---------------------------------------------------------------------------
# 3) Entropically-gated chain  (Sec. 5-6)
#    m gated stages; each: X_i -> X_{i+1} (rate p), X_i <-> Y_i (f, b).
#    Fast gating -> CV^2 -> 1/m independent of compensation strength.
# ---------------------------------------------------------------------------
def gated_chain_subgen(T, m, p0=1.0, E_p=1.0, B=800.0, K1=40.0, dH=12.0):
    """Sub-generator of an m-stage entropically-gated chain.

    K(T) = K1 * exp(dH*(1 - 1/T)) so that K(1) = K1 (trap depth at T=1) and K
    INCREASES with T (equivalently K = K0 exp(-dH/T), K0 = K1 exp(dH)); the trap
    becomes more populated as T rises, throttling the effective forward rate ->
    period-lengthening. B is the (fast) trap-return prefactor; larger B -> smaller eps.
    """
    ns = 2 * m
    Q = np.zeros((ns, ns))
    p = arrhenius(p0, E_p, T)
    b = B                               # fast, T-independent return
    K = K1 * np.exp(dH * (1.0 - 1.0 / T))
    f = b * K
    for i in range(m):
        X, Y = 2 * i, 2 * i + 1
        if i < m - 1:
            Q[X, 2 * (i + 1)] += p
        Q[X, X] -= p                    # last stage advances to absorbing
        Q[X, Y] += f
        Q[X, X] -= f
        Q[Y, X] += b
        Q[Y, Y] -= b
    return Q


def gated_chain_cv2_formula(m, eps, K):
    """Exact finite-equilibration law (Eq. 13): CV^2 = (1/m)(1 + 2 eps K/(1+K))."""
    return (1.0 / m) * (1.0 + 2.0 * eps * K / (1.0 + K))


def activity_identity_rhs(K):
    """RHS of the exact activity identity (Eq. 14): 4 K^2 / (1+K)^2 (<= 4)."""
    return 4.0 * K ** 2 / (1.0 + K) ** 2


# ---------------------------------------------------------------------------
# 4) One-trap construction for the sandwich upper bound  (Sec. 7, Fig. 5)
#    (n-1)-stage chain + a single entropic trap on stage 0, n states total.
#    Equal stage means -> CV^2 -> 1/(n-1) as eps -> 0.
# ---------------------------------------------------------------------------
def one_trap_construction_subgen(T, n, p1=1.0, E_p=0.5, B=2000.0,
                                 K1=40.0, dH=12.0):
    """Sub-generator: states X_0..X_{n-2} (chain) + Y (trap on X_0), n total."""
    ns = n
    Q = np.zeros((ns, ns))
    p1v = arrhenius(p1, E_p, T)
    b = B
    K = K1 * np.exp(dH * (1.0 - 1.0 / T))   # increases with T (period-lengthening)
    f = b * K
    cplain = p1v / (1.0 + K)            # equal-mean plain stages
    # chain X_0 -> X_1 -> ... -> X_{n-2} -> exit
    Q[0, 1] += p1v
    Q[0, 0] -= p1v
    for i in range(1, n - 2):
        Q[i, i + 1] += cplain
        Q[i, i] -= cplain
    Q[n - 2, n - 2] -= cplain            # last plain stage -> absorbing
    # entropic trap on X_0
    Q[0, n - 1] += f
    Q[0, 0] -= f
    Q[n - 1, 0] += b
    Q[n - 1, n - 1] -= b
    return Q
