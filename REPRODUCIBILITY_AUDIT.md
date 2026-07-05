# Reproducibility audit

Audit performed in the ChatGPT sandbox after cleaning and patching the repository.

## Commands run

```bash
make test
make optimum
make figures
make paper
make all
```

## Results

- `make test`: 12 tests passed.
- `make optimum`: verified the cached Figure 5 candidate in `data/optimum_n3_R2.json`.
  - `CV^2 = 0.42592661987083447`
  - analytic `R = 2.0000000000298943`
  - finite-difference audits:
    - `dT=1e-2`: `R = 2.0026001233752724`
    - `dT=3e-3`: `R = 2.0002340016302655`
    - `dT=1e-3`: `R = 2.0000260001150885`
  - sandwich bounds for `n=3`: `[1/3, 1/2]`
- `make figures`: regenerated Figures 1--5.
- `make paper`: compiled `paper/temperature_compensation.pdf` successfully with no unresolved references or LaTeX warnings in the final pass.
- `make all`: now runs tests, the cached optimum audit, figure generation, and manuscript compilation.

## Changes made during cleaning

- Removed `.pytest_cache`, `__pycache__`, and LaTeX build artifacts.
- Added `src/tempcomp/__init__.py`.
- Added `src/tempcomp/optimum.py` for deterministic verification of the cached Figure 5 candidate.
- Added `data/optimum_n3_R2.json` with the audited numerical candidate.
- Added `scripts/verify_optimum.py` and `tests/test_optimum_cache.py`.
- Updated `scripts/fig5_sandwich.py` so the Figure 5 star is loaded from the audited JSON file rather than hardcoded.
- Updated `Makefile`: `make optimum` is now a fast verification target; the slow exploratory search is `make search-optimum`.
- Updated manuscript wording to avoid overclaiming the numerical point as a certified global optimum.
- Added `CITATION.cff` and synchronized README/result numbering with the revised manuscript.


## Final manuscript audit

- In-text citation keys: 20 unique keys; all have matching `\bibitem` entries.
- Reference list entries: 20; all are cited in the manuscript.
- Cross-references: all `\ref` and `\eqref` targets resolve; no duplicate labels.
- Figures: Figures 1--5 were regenerated from scripts and visually checked in the compiled PDF render.
- Fig. 5 numerical point: loaded from `data/optimum_n3_R2.json` and audited deterministically; it is described as an upper-bound candidate, not a certified global optimum.
- Repository DOI is intentionally left as `[to be added]` in the manuscript for manual insertion.
