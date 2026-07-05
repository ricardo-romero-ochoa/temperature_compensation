"""Result 3: the entropically-gated chain achieves CV^2 = 1/m essentially
independent of the required compensation strength R (fast-gating limit)."""
import numpy as np
from tempcomp import fpt
from tempcomp.oscillators import gated_chain_subgen
from tempcomp.analysis import dlnmu_dT


def test_cv2_approaches_one_over_m():
    for m in [1, 2, 4, 6]:
        Q = gated_chain_subgen(1.0, m, p0=1.0, E_p=1.0, B=2000.0, K1=40.0, dH=12.0)
        assert abs(fpt.cv2(Q) - 1.0 / m) < 5e-3


def test_cv2_independent_of_R():
    """At fixed m, sweeping the trap enthalpy dH changes R substantially while
    CV^2 stays pinned near 1/m -- sensitivity and coherence are decoupled."""
    m = 4
    Rs, cvs = [], []
    for dH in [2.0, 6.0, 10.0]:
        subgen = lambda T: gated_chain_subgen(T, m, p0=1.0, E_p=0.5,
                                              B=3000.0, K1=60.0, dH=dH)
        Rs.append(dlnmu_dT(subgen, 1.0))
        cvs.append(fpt.cv2(subgen(1.0)))
    # R varies over a wide range ...
    assert max(Rs) - min(Rs) > 1.0
    # ... but every CV^2 sits within 5% of 1/m
    for cv in cvs:
        assert abs(cv - 1.0 / m) < 0.05 / m + 5e-3
