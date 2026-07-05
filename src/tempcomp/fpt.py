"""
fpt.py -- first-passage-time (phase-type) statistics.

A phase-type dwell is the absorption time of a CTMC with transient states and
one absorbing exit. It is specified by the transient sub-generator Q (n x n):
    Q[i, j] = rate (i -> j) for i != j (transient->transient),
    Q[i, i] = -(sum of all out-rates of i, INCLUDING the rate to the
               absorbing exit).
The exit rates are implicit in the negative row sums of Q.

Moments of the absorption time from entry state `entry` (row convention):
    m1 = ((-Q)^{-1} 1)[entry]
    m2 = 2 ((-Q)^{-2} 1)[entry]
See Neuts (1981); Aldous & Shepp (1987).
"""
from __future__ import annotations
import numpy as np


def fpt_moments(Q, entry=0):
    """Return (mean, second_moment) of the absorption time from `entry`."""
    n = Q.shape[0]
    one = np.ones(n)
    N = np.linalg.inv(-Q)                 # fundamental matrix
    m1 = (N @ one)[entry]
    m2 = 2.0 * (N @ N @ one)[entry]
    return m1, m2


def cv2(Q, entry=0):
    """Squared coefficient of variation CV^2 = Var / mean^2 of the dwell."""
    m1, m2 = fpt_moments(Q, entry)
    var = m2 - m1 * m1
    return var / (m1 * m1)


def coherence_from_cv2(cv2_value):
    """Renewal (phase-diffusion) coherence N = 1 / (pi CV^2)."""
    return 1.0 / (np.pi * cv2_value)


def erlang_subgen(n, rate=1.0):
    """Sub-generator of an order-n Erlang dwell (n identical sequential stages).
    Aldous-Shepp minimizer; CV^2 = 1/n exactly."""
    Q = np.zeros((n, n))
    for i in range(n):
        if i < n - 1:
            Q[i, i + 1] += rate
        Q[i, i] -= rate               # last stage's rate leaks to absorbing
    return Q


def sequential_subgen(rates):
    """Sub-generator of a sequence of exponential stages with given forward
    rates (hypoexponential dwell)."""
    n = len(rates)
    Q = np.zeros((n, n))
    for i, r in enumerate(rates):
        if i < n - 1:
            Q[i, i + 1] += r
        Q[i, i] -= r
    return Q


def trap_step_subgen(p, f, b):
    """Sub-generator of ONE entropically-gated step: transient states X, Y.
        X -> exit (advance)  : rate p
        X -> Y (enter trap)  : rate f
        Y -> X (leave trap)  : rate b
    Exact single-step result: CV^2 = 1 + 2 f p / (b + f)^2 = 1 + 2 eps K/(1+K),
    with K = f/b and eps = p/(f+b).
    """
    # states: 0 = X (entry), 1 = Y
    return np.array([[-(p + f), f],
                     [b, -b]], dtype=float)


def trap_step_activity(p, f, b):
    """Expected number of trap traversals (X<->Y) per step completion = 2 f / p."""
    return 2.0 * f / p
