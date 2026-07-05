"""Audit the cached Figure 5 numerical candidate."""
from pathlib import Path

from tempcomp.optimum import audit_candidate, load_candidate


def test_cached_optimum_candidate():
    path = Path(__file__).resolve().parents[1] / "data" / "optimum_n3_R2.json"
    audit = audit_candidate(load_candidate(path))
    assert 1.0 / 3.0 <= audit["CV2"] <= 1.0 / 2.0
    assert abs(audit["CV2"] - 0.42592661987083447) < 1e-12
    assert abs(audit["R_analytic"] - 2.0) < 1e-6
    assert abs(audit["R_fd_1e-3"] - 2.0) < 1e-3
