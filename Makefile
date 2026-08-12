# rietveld-agent -- reproducible Rietveld QPA (SRM 2686a)
# Professional research-workflow entry points.  See README.md and
# docs/rqpa_protocol.md for details.

PY     := .venv/bin/python
TECTONIC ?= tectonic

.PHONY: env paper figures report check cite clean

## env        -- create the local virtualenv with the runtime dependencies
env:
	python3 -m venv .venv
	$(PY) -m pip install --upgrade pip
	$(PY) -m pip install numpy scipy matplotlib

## paper      -- compile the manuscript (paper/main.tex -> paper/main.pdf)
paper:
	cd paper && $(TECTONIC) main.tex

## figures    -- regenerate paper figures from the spike-15 results
figures: report
	$(PY) paper/figures/make_figures.py

## report     -- rerun the full RQPA protocol suite (GSAS-II, ~45 min)
report:
	$(PY) benchmarks/spikes/spike_15_rqpa_protocol.py

## check      -- syntax-check every benchmark/eval module
check:
	$(PY) -m py_compile benchmarks/eval/*.py benchmarks/spikes/spike_1[4-5]_*.py
	@echo "check OK"

## cite       -- show the recommended citation (BibTeX)
cite:
	@cat CITATION.cff | sed -n '/preferred-citation/,/^[a-z]/p'

## clean      -- remove build artifacts
clean:
	rm -f paper/main.pdf paper/main.aux paper/main.bbl paper/main.blg \
	      paper/main.log paper/main.out
