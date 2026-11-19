.PHONY: help install test benchmark egfr lint clean

help:
	@grep -E "^[a-z-]+:.*##" $(MAKEFILE_LIST) | sed "s/:.*## / - /"

install:            ## editable install into the active environment
	pip install -e ".[dev]"

test:               ## unit tests (fast)
	python -m pytest -q -m "not slow"

benchmark:          ## reproduce the EGFR case study and check 16 published values
	cd examples/egfr && evodyn benchmark

egfr:               ## run the full pipeline on the shipped EGFR configuration
	cd examples/egfr && evodyn run -c config.yaml

all-tests:          ## unit tests and the EGFR regression
	python -m pytest -q

clean:
	rm -rf examples/*/run_* .pytest_cache **/__pycache__
