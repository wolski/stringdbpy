VENV_BIN := .venv/bin
WHEEL = $(firstword $(wildcard dist/*.whl))
SDIST = $(firstword $(wildcard dist/*.tar.gz))

DOCKER_IMAGE_LOCAL := string-gsea:local
DOCKER_IMAGE_REMOTE := ghcr.io/wolski/string-gsea:latest

.DEFAULT_GOAL := help
.PHONY: help sync format format-check lint typecheck deps test build check clean \
	test-smoke test-integration docker-build docker-build-local test-docker \
	render-docker render clean-integration

help:  ## Show developer commands
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

sync:  ## Synchronize the locked development environment
	uv sync --frozen --group dev

format:  ## Format and autofix source and tests
	$(VENV_BIN)/ruff format src tests
	$(VENV_BIN)/ruff check --fix src tests

format-check:  ## Check formatting without changing files
	$(VENV_BIN)/ruff format --check src tests

lint:  ## Run code and import-architecture lint checks
	$(VENV_BIN)/ruff check src tests
	$(VENV_BIN)/lint-imports

typecheck:  ## Run standard Pyright in strict mode
	$(VENV_BIN)/pyright

deps:  ## Validate dependency declarations
	$(VENV_BIN)/deptry .

test:  ## Run unit tests with branch coverage
	$(VENV_BIN)/pytest --cov --cov-branch

build:  ## Build and validate source and wheel distributions
	uv build --clear
	$(VENV_BIN)/twine check dist/*
	wheel_uri=$$($(VENV_BIN)/python -c \
		'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve().as_uri())' "$(WHEEL)"); \
	wheel_hash=$$($(VENV_BIN)/python -c \
		'import hashlib, sys; from pathlib import Path; print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())' "$(WHEEL)"); \
	uv run --isolated --no-project --with "string-gsea @ $$wheel_uri#sha256=$$wheel_hash" python -c \
		'from importlib import metadata, resources; import shutil; \
		assert metadata.version("string-gsea"); \
		assert resources.files("string_gsea.workflow").joinpath("Snakefile").is_file(); \
		assert resources.files("string_gsea.data.mappings").joinpath("NCBI_nodes.zip").is_file(); \
		assert resources.files("string_gsea.ora").joinpath("enrichment_results.schema.json").is_file(); \
		assert all(shutil.which(command) for command in ("string_gsea_run", "string_ora_run", "string_gsea_workflow", "_string_gsea_render_only"))'
	sdist_uri=$$($(VENV_BIN)/python -c \
		'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve().as_uri())' "$(SDIST)"); \
	sdist_hash=$$($(VENV_BIN)/python -c \
		'import hashlib, sys; from pathlib import Path; print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())' "$(SDIST)"); \
	uv run --isolated --no-project --with "string-gsea @ $$sdist_uri#sha256=$$sdist_hash" python -c \
		'from importlib import metadata, resources; \
		assert metadata.version("string-gsea"); \
		assert resources.files("string_gsea.workflow").joinpath("Snakefile").is_file(); \
		assert resources.files("string_gsea.data.mappings").joinpath("species.v12.0.zip").is_file(); \
		assert resources.files("string_gsea.ora").joinpath("enrichment_results.schema.json").is_file()'
	wheel_uri=$$($(VENV_BIN)/python -c \
		'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve().as_uri())' "$(WHEEL)"); \
	wheel_hash=$$($(VENV_BIN)/python -c \
		'import hashlib, sys; from pathlib import Path; print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())' "$(WHEEL)"); \
	uv run --isolated --no-project --with "string-gsea @ $$wheel_uri#sha256=$$wheel_hash" sh -c \
		'for command in string_gsea_write_config string_gsea_run string_ora_run string_gsea_workflow \
		_string_gsea_render _string_gsea_render_only _string_gsea_package; do "$$command" --help >/dev/null; done'

check:  ## Run every merge-blocking quality gate
	uv lock --check
	$(MAKE) format-check lint typecheck deps test build

clean:  ## Remove generated build and quality artifacts
	$(VENV_BIN)/python -c "import shutil; [shutil.rmtree(path, ignore_errors=True) for path in ('build', 'dist', '.pytest_cache', '.ruff_cache')]"

test-smoke:  ## Run the quick workflow test (STRING-DB + Quarto + R)
	$(VENV_BIN)/pytest -m smoke tests -v -s

test-integration:  ## Run the full workflow test suite (STRING-DB + Quarto + R)
	$(VENV_BIN)/pytest -m integration tests -v -s

docker-build-local:  ## Build the local amd64 Docker image used by CI
	docker buildx build --platform linux/amd64 -f docker/Dockerfile -t $(DOCKER_IMAGE_LOCAL) --load .

docker-build:  ## Build the GHCR-tagged Docker image (requires a GHCR login)
	docker buildx build -f docker/Dockerfile -t $(DOCKER_IMAGE_REMOTE) --load .

test-docker:  ## Run the mouse XLSX workflow in Docker
	rm -rf tests/data/outputs/mouse_xlsx_docker
	./docker/string_gsea_docker.sh --image-repo string-gsea --image-version local \
		tests/data/datasets/mouse_xlsx/input.zip mouse_fasta \
		tests/data/outputs/mouse_xlsx_docker \
		--which pep_2_no_imputed --cores 1
	@test -f tests/data/outputs/mouse_xlsx_docker/outputs.yml
	@test -f tests/data/outputs/mouse_xlsx_docker/WU_mouse_fasta_GSEA.zip

render-docker:  ## Re-render reports for existing Docker results
	./docker/string_gsea_docker.sh --image-repo string-gsea --image-version local render \
		tests/data/outputs/mouse_xlsx_docker mouse_fasta

render:  ## Re-render reports locally
	uv run --frozen _string_gsea_render_only tests/data/outputs/mouse_xlsx_docker mouse_fasta

clean-integration:  ## Remove generated integration-test outputs
	rm -rf tests/data/outputs/
