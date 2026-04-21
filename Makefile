DOCKER_IMAGE_LOCAL := string-gsea:local
DOCKER_IMAGE_REMOTE := ghcr.io/wolski/string-gsea:latest

.DEFAULT_GOAL := help
.PHONY: help check test-smoke test-integration docker-build docker-build-local test-docker render-docker clean-integration lint format

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

docker-build-local:            ## Build Docker image locally as string-gsea:local (amd64, matches CI)
	docker buildx build --platform linux/amd64 -f docker/Dockerfile -t $(DOCKER_IMAGE_LOCAL) --load .

docker-build:                  ## Build Docker image as ghcr.io tag (CI only — requires docker login ghcr.io)
	docker buildx build -f docker/Dockerfile -t $(DOCKER_IMAGE_REMOTE) --load .

test-docker:                   ## Run mouse_xlsx workflow in Docker (requires built image)
	rm -rf tests/data/outputs/mouse_xlsx_docker
	./docker/string_gsea_docker.sh --image-repo string-gsea --image-version local \
		tests/data/datasets/mouse_xlsx/input.zip mouse_fasta \
		tests/data/outputs/mouse_xlsx_docker \
		--which pep_2_no_imputed --cores 1
	@echo "--- Verifying outputs ---"
	@test -f tests/data/outputs/mouse_xlsx_docker/outputs.yml && echo "PASS: outputs.yml exists" || (echo "FAIL: outputs.yml missing" && exit 1)
	@test -f tests/data/outputs/mouse_xlsx_docker/WU_mouse_fasta_GSEA.zip && echo "PASS: zip exists" || (echo "FAIL: zip missing" && exit 1)

render-docker:                 ## Re-render reports for existing mouse_xlsx_docker results (no STRING re-run)
	./docker/string_gsea_docker.sh --image-repo string-gsea --image-version local render \
		tests/data/outputs/mouse_xlsx_docker mouse_fasta

render:                        ## Re-render reports locally (no Docker, no STRING re-run)
	uv run _string_gsea_render_only tests/data/outputs/mouse_xlsx_docker mouse_fasta

clean-integration:             ## Remove integration test outputs (forces re-run)
	rm -rf tests/data/outputs/
