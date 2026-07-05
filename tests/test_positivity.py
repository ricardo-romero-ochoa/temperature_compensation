"""Proposition 1 (positivity lemma): a purely sequential Arrhenius clock cannot
be temperature compensated -- period strictly decreases with T (Q10 > 1) and all
control coefficients are negative."""
import numpy as np
from tempcomp.analysis import dlntau_dT_sequential, control_coefficients_sequential


def _rates(T, As, Es):
    return [A * np.exp(-E / T) for A, E in zip(As, Es)]


def test_sequential_cannot_compensate():
    As = [1.0, 0.7, 1.3, 0.9]
    Es = [4.0, 8.0, 2.0, 11.0]      # arbitrary positive activation energies
    rates_fn = lambda T: _rates(T, As, Es)
    # d ln tau / dT < 0 strictly (period speeds up with T -> no compensation)
    assert dlntau_dT_sequential(rates_fn, 1.0) < -1e-6


def test_all_control_coefficients_negative():
    As = [1.0, 0.5, 2.0]
    Es = [3.0, 9.0, 6.0]
    rates_fn = lambda T: _rates(T, As, Es)
    C = control_coefficients_sequential(rates_fn, 1.0)
    assert np.all(C < 0)
    # summation theorem: control coefficients of the period sum to -1
    assert abs(C.sum() + 1.0) < 1e-6
