"""
analysis.py -- derived quantities: compensation deficit, temperature
sensitivity of first-passage means, and finite-difference helpers.
"""
from __future__ import annotations
import numpy as np
from . import markov, fpt


def omega_of_T_two_loop(edges_fn, T, **kw):
    """omega(T) = |Im lambda| of the two-loop model's clock mode."""
    from .markov import generator, eigen_coherence
    from .oscillators import N_STATES_TWO_LOOP
    L = generator(edges_fn(T, **kw), N_STATES_TWO_LOOP)
    return eigen_coherence(L)[0]


def compensation_deficit_omega(edges_fn, T0, dT=1e-3, **kw):
    """chi = |d ln omega / dT| at T0 (central difference), eigenvalue omega."""
    op = omega_of_T_two_loop(edges_fn, T0 + dT, **kw)
    om = omega_of_T_two_loop(edges_fn, T0 - dT, **kw)
    if op <= 0 or om <= 0:
        return np.nan
    return abs((np.log(op) - np.log(om)) / (2 * dT))


def dlnmu_dT(subgen_fn, T0, entry=0, dT=1e-3):
    """R = d ln(mean FPT) / dT for a temperature-dependent sub-generator
    builder subgen_fn(T) -> Q. Central difference."""
    mp, _ = fpt.fpt_moments(subgen_fn(T0 + dT), entry)
    mm, _ = fpt.fpt_moments(subgen_fn(T0 - dT), entry)
    m0, _ = fpt.fpt_moments(subgen_fn(T0), entry)
    return (mp - mm) / (2 * dT) / m0


def dlntau_dT_sequential(rates_fn, T0, dT=1e-3):
    """d ln(tau)/dT for a sequential clock tau(T) = sum 1/k_i(T),
    rates_fn(T) -> list of forward rates. Used for the positivity lemma."""
    def tau(T):
        return sum(1.0 / k for k in rates_fn(T))
    return (np.log(tau(T0 + dT)) - np.log(tau(T0 - dT))) / (2 * dT)


def control_coefficients_sequential(rates_fn, T0, dlnk=1e-4):
    """C_j = d ln tau / d ln k_j for a sequential clock at T0.
    Returns array of control coefficients (all negative for a sequential
    clock -- the positivity lemma)."""
    ks = np.array(rates_fn(T0), dtype=float)
    def tau_from(ks_):
        return sum(1.0 / k for k in ks_)
    tau0 = tau_from(ks)
    C = np.zeros(len(ks))
    for j in range(len(ks)):
        kp = ks.copy(); kp[j] *= np.exp(dlnk)
        km = ks.copy(); km[j] *= np.exp(-dlnk)
        C[j] = (np.log(tau_from(kp)) - np.log(tau_from(km))) / (2 * dlnk)
    return C
