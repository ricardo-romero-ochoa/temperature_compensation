# Reproduction pipeline for the temperature-compensation results.
PY ?= python3

.PHONY: help install test figures optimum search-optimum paper all clean

help:
	@echo "targets:"
	@echo "  make install        - install the package (editable) + deps"
	@echo "  make test           - run the pytest suite (asserts exact analytic claims)"
	@echo "  make figures        - regenerate Figures 1-5 into figures/"
	@echo "  make optimum        - fast audit of the cached Fig. 5 numerical candidate"
	@echo "  make search-optimum - stochastic constrained optimization search (slow)"
	@echo "  make paper          - compile the LaTeX paper (needs pdflatex)"
	@echo "  make all            - test + optimum audit + figures + paper"
	@echo "  make clean          - remove generated figures and LaTeX build artifacts"

install:
	$(PY) -m pip install -e ".[test]"

test:
	$(PY) -m pytest tests/ -q

figures:
	$(PY) scripts/fig1_twoloop_tradeoff.py
	$(PY) scripts/fig2_closedform_floor.py
	$(PY) scripts/fig3_nofloor.py
	$(PY) scripts/fig4_finite_eps.py
	$(PY) scripts/fig5_sandwich.py

optimum:
	$(PY) scripts/verify_optimum.py

search-optimum:
	$(PY) scripts/optimality_search.py 3 2.0

paper: figures
	cd paper && pdflatex -interaction=nonstopmode temperature_compensation.tex && \
	  pdflatex -interaction=nonstopmode temperature_compensation.tex

all: test optimum figures paper

clean:
	rm -f figures/*.png figures/*.pdf figures/*.svg
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.bbl paper/*.blg
