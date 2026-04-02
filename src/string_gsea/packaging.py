"""Helpers for packaging STRING-GSEA results for delivery."""

from __future__ import annotations

from pathlib import Path

import yaml


def write_outputs_yml(output_file: str, zip_path: str) -> None:
    """Write outputs.yml for bfabric-app-runner staging.

    Follows the same pattern as diann_runner's write_outputs_yml.
    """
    outputs = [
        {
            "local_path": str(Path(zip_path).resolve()),
            "store_entry_path": Path(zip_path).name,
            "type": "bfabric_copy_resource",
        },
    ]
    data = {"outputs": outputs}
    with open(output_file, "w") as f:
        yaml.dump(data, f, default_flow_style=False)
