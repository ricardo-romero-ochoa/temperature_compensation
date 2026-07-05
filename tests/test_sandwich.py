"""Theorem 1 (sandwich): 1/n <= CV^2_min(n,R) <= 1/(n-1).
Lower bound realized exactly by Erlang (Aldous-Shepp); upper bound realized by
the one-trap construction."""
import numpy as np
from tempcomp import fpt
from tempcomp.oscillators import one_trap_construction_subgen


def test_erlang_is_one_over_n():
    """Aldous-Shepp minimizer: Erlang-n has CV^2 = 1/n exactly."""
    for n in [2, 3, 4, 5, 8]:
        assert abs(fpt.cv2(fpt.erlang_subgen(n, 1.3)) - 1.0 / n) < 1e-10


def test_construction_in_sandwich():
    """The compensating one-trap construction lies in [1/n, 1/(n-1)]."""
    for n in [3, 4, 5, 6]:
        cv = fpt.cv2(one_trap_construction_subgen(1.0, n))
        assert 1.0 / n - 1e-6 <= cv <= 1.0 / (n - 1) + 1e-3


def test_construction_beats_or_meets_upper_bound():
    """Construction achieves (at most) the 1/(n-1) upper bound."""
    for n in [3, 4, 5, 6]:
        cv = fpt.cv2(one_trap_construction_subgen(1.0, n))
        assert cv <= 1.0 / (n - 1) + 1e-3
