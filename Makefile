.DEFAULT_GOAL := help
.PHONY: help test test-integration test-all test-ci clean-ci lint format check

help:                          ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

## Testing
test:                          ## Run fast/unit tests only (default)
	uv run pytest tests

test-integration:              ## Run integration tests (requires STRING-DB access)
	uv run pytest -m integration tests

test-all:                      ## Run all tests (unit + integration)
	uv run pytest -m "" tests

test-ci:                       ## Run Snakemake CI datasets (small, single-contrast)
	cd tests && uv run snakemake -s Snakefile ci -j1

clean-ci:                      ## Remove Snakemake CI outputs (forces re-run)
	cd tests && uv run snakemake -s Snakefile clean -j1

## Code quality
lint:                          ## Check linting
	uv run ruff check src/ tests/

format:                        ## Auto-format code
	uv run ruff format src/ tests/

check: lint test               ## Lint + unit tests
