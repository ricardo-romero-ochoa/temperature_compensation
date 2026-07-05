[![DOI](https://zenodo.org/badge/1289677055.svg)](https://doi.org/10.5281/zenodo.21200894)


# tempcomp

Reproducibility package for:

> **The precision cost of temperature compensation in biochemical clocks: internal structure and dynamical activity, not entropy production**  
> Ricardo Romero

This repository is designed so that the analytical identities, numerical figures,
and the cached Figure 5 optimization candidate can be regenerated or audited from
source. The rigorous results of the paper do **not** depend on the numerical
candidate; it is included as an audited, reproducible point inside the sandwich.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```

Tested with Python 3.12. The package requires `numpy`, `scipy`, `sympy`,
`matplotlib`, and `pytest`.

## Reproduce the package

```bash
make test           # pytest suite for the exact analytical claims
make optimum        # fast audit of data/optimum_n3_R2.json
make figures        # regenerate Figures 1-5 into figures/
make all            # test + optimum audit + figures + manuscript PDF
make paper          # regenerate figures and compile the manuscript PDF; requires pdflatex
```

The slow stochastic optimizer is kept separate:

```bash
make search-optimum # exploratory search for the n=3, R=2 candidate
```

## Repository layout

```text
src/tempcomp/
  __init__.py       package marker and version
  markov.py         finite CTMC generators, steady states, eigenvalue coherence,
                   entropy production, dynamical activity
  fpt.py            phase-type first-passage moments, CV^2, Erlang and gated steps
  oscillators.py    two-loop model, semi-Markov closed forms, gated chain,
                   one-trap sandwich construction
  analysis.py       compensation deficit and temperature-sensitivity helpers
  optimum.py        cached-candidate loading and analytic/finite-difference audit
scripts/            one script per figure + optimum verification/search
src/                importable Python package
figures/            generated figures
paper/              LaTeX source and compiled PDF
data/               cached audited Fig. 5 numerical candidate
tests/              pytest suite
```

## Result-to-code map

| Paper item | Reproduced or audited by |
|---|---|
| Proposition 1: sequential Arrhenius clocks cannot compensate | `tests/test_positivity.py` |
| Fig. 1: canonical two-loop trade-off; `chi -> 0`, coherence falls, entropy production does not rise | `scripts/fig1_twoloop_tradeoff.py` |
| Eqs. 5-8 / Fig. 2: semi-Markov closed forms, `1/pi` floor | `tests/test_semimarkov.py`, `scripts/fig2_closedform_floor.py` |
| Result 2 / Fig. 3: gated-chain `CV^2 -> 1/m` independent of compensation strength | `tests/test_gated_chain.py`, `scripts/fig3_nofloor.py` |
| Result 3 / Fig. 4: finite-equilibration law and activity identity | `tests/test_finite_eps.py`, `scripts/fig4_finite_eps.py` |
| Theorem 1 / Fig. 5: sandwich `1/n <= CV^2_min <= 1/(n-1)` | `tests/test_sandwich.py`, `scripts/fig5_sandwich.py` |
| Fig. 5 star: audited numerical candidate at `n=3, R=2` | `data/optimum_n3_R2.json`, `scripts/verify_optimum.py`, `tests/test_optimum_cache.py` |

## Cached optimization candidate

The file `data/optimum_n3_R2.json` stores a numerical candidate for the order-3
Arrhenius phase-type search. `make optimum` verifies that it:

- lies inside the sandwich `[1/3, 1/2]`,
- has `CV^2 = 0.42592661987083447`,
- has analytic sensitivity `R = 2.0000000000298943`, and
- passes central-difference audits at `dT = 1e-2, 3e-3, 1e-3`.

This candidate is an audited upper-bound point, not a certificate of global
optimality. The exact constrained optimum remains the open problem stated in the
paper.

## Conventions

- Units: `k_B = 1`.
- Arrhenius rates: `k(T) = k(T0) exp(E(1/T0 - 1/T))` when rates are stored at `T0`.
- Finite CTMC generator: `dp/dt = L p`, with `L[i,j] = rate(j -> i)` and columns summing to zero.
- Phase-type sub-generator `Q`: row convention; exit rates are implicit in the negative row sums.
- Renewal coherence: `N = 1/(pi CV^2)`.
- Finite two-loop coherence: selected oscillatory eigenmode, `|Im lambda|/|Re lambda|`.

## Scope and limitations

- `CV^2 -> 1/m` is a fast-gating limit; finite equilibration gives the exact correction tested in `tests/test_finite_eps.py`.
- The renewal treatment assumes independent stage dwell times.
- The sandwich is rigorous on both sides but is not tight. The cached numerical point is reproducible evidence that the upper construction is not saturated for one parameter set, not a proof of the exact optimum.

## License

MIT; see `LICENSE`.
