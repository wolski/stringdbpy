#!/bin/bash
set -euxo pipefail
task_dir="$1"
cp /home/bfabric/slurmworker/config/A366_DIANN/run.sh "$task_dir/.."
ln -s /home/bfabric/slurmworker/config/A366_DIANN/Snakefile "$task_dir/.."
cp /home/bfabric/slurmworker/config/A366_DIANN/uv.lock "$task_dir/.."
cp /home/bfabric/slurmworker/config/A366_DIANN/pyproject.toml "$task_dir/.."

cd "$task_dir/.."
uv run --verbose --locked --project "." snakemake -d "$task_dir" -p  all
