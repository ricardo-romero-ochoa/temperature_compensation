"""Semi-Markov (single-trap) reduction: exact node-0 waiting-time moments
(Eq. 5) verified symbolically, and the closed forms of Eq. (6)."""
import numpy as np
import sympy as sp
from tempcomp import fpt


def test_psi0_moments_symbolic():
    s, a, d, gam = sp.symbols("s a d gamma", positive=True)
    psi_e = gam / (gam + s)                       # exponential excursion
    psi0 = a / (a + s + d * (1 - psi_e))          # renewal-equation solution
    mean = sp.simplify(-sp.diff(psi0, s).subs(s, 0))
    m2 = sp.simplify(sp.diff(psi0, s, 2).subs(s, 0))
    var = sp.simplify(m2 - mean ** 2)
    # in units a=1, mu_e = 1/gamma, d = n_ex:  mean = 1 + n_ex mu_e,
    # var = 1 + 2 n_ex mu_e + n_ex(2+n_ex) mu_e^2
    mu_e, n_ex = sp.symbols("mu_e n_ex", positive=True)
    subs = {gam: 1 / mu_e, a: 1, d: n_ex}
    mean_u = sp.simplify(mean.subs(subs))
    var_u = sp.simplify(var.subs(subs))
    assert sp.simplify(mean_u - (1 + n_ex * mu_e)) == 0
    assert sp.simplify(var_u - (1 + 2 * n_ex * mu_e
                                + n_ex * (2 + n_ex) * mu_e ** 2)) == 0


def test_single_trap_cv2_numeric():
    # single gated step: CV^2 = 1 + 2 eps K/(1+K), K=f/b, eps=p/(f+b)
    for p, f, b in [(1.0, 20.0, 1.0), (0.3, 60.0, 3.0), (2.0, 5.0, 4.0)]:
        K, eps = f / b, p / (f + b)
        cv = fpt.cv2(fpt.trap_step_subgen(p, f, b))
        assert abs(cv - (1 + 2 * eps * K / (1 + K))) < 1e-10
