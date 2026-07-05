"""Utilities for auditing cached constrained-optimization candidates.

The theorem's sandwich bounds are analytic. The Figure 5 star is a numerical
candidate inside the sandwich; these functions verify that a cached candidate
actually satisfies the stated CV^2 and temperature-sensitivity audits.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from . import fpt


def load_candidate(path: str | Path) -> dict[str, Any]:
    """Load a cached optimization candidate from JSON."""
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def rates_at_T(candidate: dict[str, Any], T: float) -> np.ndarray:
    """Return rates at temperature T from rates stored at T0.

    Arrhenius convention: k(T)=k(T0)*exp(E*(1/T0 - 1/T)), with k_B=1.
    """
    k0 = np.asarray(candidate["rates_at_T0"], dtype=float)
    E = np.asarray(candidate["activation_energies"], dtype=float)
    T0 = float(candidate.get("T0", 1.0))
    return k0 * np.exp(E * (1.0 / T0 - 1.0 / T))


def subgenerator(candidate: dict[str, Any], T: float | None = None) -> np.ndarray:
    """Build the transient sub-generator Q(T) of a cached candidate."""
    n = int(candidate["n"])
    T = float(candidate.get("T0", 1.0) if T is None else T)
    rates = rates_at_T(candidate, T)
    edges = [tuple(x) for x in candidate["transient_edge_order"]]
    exits = list(candidate["exit_state_order"])
    if len(rates) != len(edges) + len(exits):
        raise ValueError("rates_at_T0 length does not match edge/exit order")
    Q = np.zeros((n, n), dtype=float)
    for idx, (i, j) in enumerate(edges):
        r = rates[idx]
        Q[i, j] += r
        Q[i, i] -= r
    offset = len(edges)
    for j, i in enumerate(exits):
        Q[i, i] -= rates[offset + j]
    return Q


def subgenerator_derivative(candidate: dict[str, Any], T: float | None = None) -> np.ndarray:
    """Analytic dQ/dT under the Arrhenius convention used by ``rates_at_T``."""
    n = int(candidate["n"])
    T = float(candidate.get("T0", 1.0) if T is None else T)
    rates = rates_at_T(candidate, T)
    E = np.asarray(candidate["activation_energies"], dtype=float)
    dr = rates * E / (T * T)
    edges = [tuple(x) for x in candidate["transient_edge_order"]]
    exits = list(candidate["exit_state_order"])
    Qp = np.zeros((n, n), dtype=float)
    for idx, (i, j) in enumerate(edges):
        r = dr[idx]
        Qp[i, j] += r
        Qp[i, i] -= r
    offset = len(edges)
    for j, i in enumerate(exits):
        Qp[i, i] -= dr[offset + j]
    return Qp


def cv2(candidate: dict[str, Any], T: float | None = None) -> float:
    """Squared coefficient of variation of the candidate dwell."""
    return float(fpt.cv2(subgenerator(candidate, T)))


def sensitivity_analytic(candidate: dict[str, Any], T: float | None = None) -> float:
    """Analytic R=d ln(mean FPT)/dT at T."""
    T = float(candidate.get("T0", 1.0) if T is None else T)
    Q = subgenerator(candidate, T)
    Qp = subgenerator_derivative(candidate, T)
    F = np.linalg.inv(-Q)
    one = np.ones(int(candidate["n"]))
    entry = int(candidate.get("state_entry", 0))
    alpha = np.zeros(int(candidate["n"]))
    alpha[entry] = 1.0
    mean = alpha @ F @ one
    dmean = alpha @ F @ Qp @ F @ one
    return float(dmean / mean)


def sensitivity_central(candidate: dict[str, Any], dT: float, T: float | None = None) -> float:
    """Central-difference audit of R=d ln(mean FPT)/dT."""
    T = float(candidate.get("T0", 1.0) if T is None else T)
    entry = int(candidate.get("state_entry", 0))
    m0, _ = fpt.fpt_moments(subgenerator(candidate, T), entry)
    mp, _ = fpt.fpt_moments(subgenerator(candidate, T + dT), entry)
    mm, _ = fpt.fpt_moments(subgenerator(candidate, T - dT), entry)
    return float((mp - mm) / (2.0 * dT) / m0)


def audit_candidate(candidate: dict[str, Any], r_tol: float = 0.01) -> dict[str, float]:
    """Return a reproducibility audit and raise ValueError if it fails.

    ``r_tol`` is loose enough for the deliberately coarse dT=1e-2 finite
    difference but strict for the analytic value and the smaller step sizes.
    """
    n = int(candidate["n"])
    target = float(candidate["R_target"])
    Q = subgenerator(candidate)
    cv = cv2(candidate)
    R = sensitivity_analytic(candidate)
    cond = float(np.linalg.cond(Q))
    lower, upper = 1.0 / n, 1.0 / (n - 1)
    diffs = {d: sensitivity_central(candidate, d) for d in (1e-2, 3e-3, 1e-3)}
    if not (lower <= cv <= upper):
        raise ValueError(f"CV^2={cv} is outside sandwich [{lower}, {upper}]")
    if abs(R - target) > 5e-6:
        raise ValueError(f"analytic R={R} does not match target {target}")
    if max(abs(v - target) for v in diffs.values()) > r_tol:
        raise ValueError(f"finite-difference R audit failed: {diffs}")
    return {
        "CV2": cv,
        "R_analytic": R,
        "condition_number": cond,
        "R_fd_1e-2": diffs[1e-2],
        "R_fd_3e-3": diffs[3e-3],
        "R_fd_1e-3": diffs[1e-3],
        "sandwich_lower": lower,
        "sandwich_upper": upper,
    }
