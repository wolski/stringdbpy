#!/bin/bash
set -euxo pipefail
task_dir="$1"
cp /home/bfabric/slurmworker/config/A373_DIANN/run.sh "$task_dir/.."
ln -s /home/bfabric/slurmworker/config/A373_DIANN/testing_string_api.py "$task_dir/.."
cp /home/bfabric/slurmworker/config/A373_DIANN/uv.lock "$task_dir/.."
cp /home/bfabric/slurmworker/config/A373_DIANN/pyproject.toml "$task_dir/.."

cd "$task_dir/.."
uv run --verbose --locked --project "." testing_string_api.py
