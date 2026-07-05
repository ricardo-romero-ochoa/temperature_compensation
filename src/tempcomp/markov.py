"""
markov.py -- finite continuous-time Markov chain machinery.

Convention: dp/dt = L p, with L[i, j] = rate (j -> i) for i != j and
L[j, j] = -sum_{i != j} rate(j -> i). Columns of L sum to zero.

Edges are given as a list of tuples (src, dst, rate).
"""
from __future__ import annotations
import numpy as np


def generator(edges, n):
    """Build an n x n generator L from an edge list [(src, dst, rate), ...]."""
    L = np.zeros((n, n))
    for s, d, r in edges:
        L[d, s] += r
        L[s, s] -= r
    return L


def steady_state(L):
    """Stationary distribution: normalized right null vector of L."""
    w, v = np.linalg.eig(L)
    i = int(np.argmin(np.abs(w)))
    pi = np.real(v[:, i])
    return pi / pi.sum()


def leading_oscillatory_mode(L):
    """Return the complex eigenvalue with the largest |Im| (the clock mode).

    Selecting by largest imaginary part tracks the oscillation robustly by
    continuity, avoiding spurious near-zero complex spectator modes that a
    'largest real part' rule would grab.
    """
    w, _ = np.linalg.eig(L)
    cplx = w[np.abs(w.imag) > 1e-9]
    if len(cplx) == 0:
        return None
    return cplx[int(np.argmax(np.abs(cplx.imag)))]


def eigen_coherence(L):
    """(omega, N) for the leading oscillatory mode.

    omega = |Im lambda|, coherence N = |Im lambda| / |Re lambda|
    (Barato-Seifert definition). Returns (0.0, 0.0) if non-oscillatory.
    """
    lam = leading_oscillatory_mode(L)
    if lam is None or lam.real == 0:
        return 0.0, 0.0
    return abs(lam.imag), abs(lam.imag) / abs(lam.real)


def entropy_production(edges, pi):
    """Steady-state entropy production sigma = sum_edges (Jf - Jb) ln(Jf/Jb),
    in units of k_B. Each undirected edge counted once."""
    kd = {(s, d): r for s, d, r in edges}
    seen, sigma = set(), 0.0
    for s, d, r in edges:
        key = frozenset((s, d))
        if key in seen:
            continue
        seen.add(key)
        jf = kd.get((s, d), 0.0) * pi[s]
        jb = kd.get((d, s), 0.0) * pi[d]
        if jf > 1e-300 and jb > 1e-300:
            sigma += (jf - jb) * np.log(jf / jb)
    return sigma


def dynamical_activity(edges, pi):
    """Total steady-state jump rate (traffic) K = sum_{edges} pi_s k_{s->d}."""
    return sum(r * pi[s] for s, d, r in edges)
