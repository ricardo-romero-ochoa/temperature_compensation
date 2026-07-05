"""Result 4: exact finite-equilibration law CV^2 = (1/m)(1 + 2 eps K/(1+K))
and the exact activity identity (CV^2_step - 1) * K_step = 4 K^2/(1+K)^2 <= 4."""
import numpy as np
from tempcomp import fpt
from tempcomp.oscillators import gated_chain_subgen, gated_chain_cv2_formula, activity_identity_rhs


def test_finite_eps_law_exact():
    """Exact 2m-state chain CV^2 equals the closed-form law to 5+ digits."""
    for m in [1, 3, 6]:
        for B in [50.0, 400.0, 3000.0]:
            for K1 in [2.0, 20.0]:
                T = 1.0
                Q = gated_chain_subgen(T, m, p0=1.0, E_p=1.0, B=B, K1=K1, dH=12.0)
                p = 1.0 * np.exp(-1.0 / T)            # p0 exp(-E_p/T)
                K = K1                                 # K(T=1)=K1 by construction
                eps = p / (B * (1 + K))                # p/(f+b), f=b K, b=B
                assert abs(fpt.cv2(Q) - gated_chain_cv2_formula(m, eps, K)) < 1e-5


def test_activity_identity_exact():
    for p, f, b in [(1.0, 20.0, 1.0), (0.02, 4.0, 2.0), (0.5, 3.0, 3.0)]:
        K = f / b
        cv = fpt.cv2(fpt.trap_step_subgen(p, f, b))
        act = fpt.trap_step_activity(p, f, b)
        lhs = (cv - 1.0) * act
        assert abs(lhs - activity_identity_rhs(K)) < 1e-8
        assert activity_identity_rhs(K) <= 4.0 + 1e-12
