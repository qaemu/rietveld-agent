# rietveld-agent -- reproducible Rietveld QPA (SRM 2686a)
# Professional research-workflow entry points.  See README.md and
# docs/rqpa_protocol.md for details.

PY     := .venv/bin/python
PIP    := .venv/bin/pip
TECTONIC ?= tectonic

.PHONY: env paper figures report check test clean

## env        -- create the local virtualenv with the runtime dependencies
env:
	python3 -m venv .venv
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

## paper      -- compile all paper PDFs (paper/*.tex -> paper/*.pdf)
paper:
	cd paper && $(TECTONIC) main.tex
	cd paper && $(TECTONIC) paper1_protocol.tex
	cd paper && $(TECTONIC) paper2_validation.tex
	cd paper && $(TECTONIC) paper3_deployment.tex

## figures    -- regenerate paper figures from the spike-15 results
figures: report
	$(PY) paper/figures/make_figures.py

## report     -- rerun the full RQPA protocol suite (GSAS-II, ~45 min)
report:
	$(PY) benchmarks/spikes/spike_15_rqpa_protocol.py

## check      -- syntax-check every python module
check:
	$(PY) -m compileall -q core cli benchmarks tests
	@echo "check OK"

## test       -- run the regression suite (pytest)
test:
	$(PY) -m pytest -q

## clean      -- remove build artifacts
clean:
	rm -f paper/main.pdf paper/main.aux paper/main.bbl paper/main.blg \
	      paper/main.log paper/main.out paper/main.toc \
	      paper/paper1_protocol.pdf paper/paper2_validation.pdf \
	      paper/paper3_deployment.pdf