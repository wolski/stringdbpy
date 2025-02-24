#!/bin/bash
set -euxo pipefail
task_dir="$1"
cp /home/bfabric/slurmworker/config/A373_STRING_GSEA/run.sh "$task_dir/.."
#ln -s /home/bfabric/slurmworker/config/A373_STRING_GSEA/testing_string_api.py "$task_dir/.."
#cp /home/bfabric/slurmworker/config/A373_STRING_GSEA/uv.lock "$task_dir/.."
#cp /home/bfabric/slurmworker/config/A373_STRING_GSEA/pyproject.toml "$task_dir/.."

uv run --verbose --locked --project /home/bfabric/slurmworker/config/A373_STRING_GSEA -m a373_string_gsea.run_string_gsea
