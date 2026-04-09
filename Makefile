DOCKER_IMAGE := ghcr.io/wolski/string-gsea:latest

.DEFAULT_GOAL := help
.PHONY: help check test-smoke test-integration docker-build test-docker clean-integration lint format

help:                          ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

## Code quality
lint:                          ## Check linting
	uv run ruff check src/ tests/

format:                        ## Auto-format code
	uv run ruff format src/ tests/

check: lint                    ## Lint + unit tests
	uv run pytest tests

## Integration
test-smoke:                    ## Quick workflow test (yeast RNK, single contrast)
	uv run pytest -m smoke tests -v -s

test-integration:              ## Run full workflow on all CI datasets (requires STRING-DB + Quarto + R)
	uv run pytest -m integration tests -v -s

docker-build:                  ## Build the Docker image locally
	docker buildx build -f docker/Dockerfile -t $(DOCKER_IMAGE) --load .

test-docker:                   ## Run mouse_xlsx workflow in Docker (requires built image)
	rm -rf tests/data/outputs/mouse_xlsx_docker
	./docker/string_gsea_docker.sh \
		tests/data/datasets/mouse_xlsx/input.zip mouse_fasta \
		tests/data/outputs/mouse_xlsx_docker \
		--which pep_2_no_imputed --cores 1
	@echo "--- Verifying outputs ---"
	@test -f tests/data/outputs/mouse_xlsx_docker/outputs.yml && echo "PASS: outputs.yml exists" || (echo "FAIL: outputs.yml missing" && exit 1)
	@test -f tests/data/outputs/mouse_xlsx_docker/WU_mouse_fasta_GSEA.zip && echo "PASS: zip exists" || (echo "FAIL: zip missing" && exit 1)

clean-integration:             ## Remove integration test outputs (forces re-run)
	rm -rf tests/data/outputs/
