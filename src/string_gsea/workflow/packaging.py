"""Packaging of completed workflow results for B-Fabric delivery."""

import shutil
from pathlib import Path

import yaml


def package_results(output_base: Path, workunit_id: str, outputs_yml_path: Path) -> Path:
    """Archive a result directory and write the established outputs manifest."""
    workunit = output_base / f"WU_{workunit_id}_GSEA"
    archive = output_base / f"WU_{workunit_id}_GSEA.zip"
    shutil.make_archive(
        str(archive).removesuffix(".zip"),
        "zip",
        root_dir=str(workunit.parent),
        base_dir=workunit.name,
    )
    document = {
        "outputs": [
            {
                "local_path": str(archive.resolve()),
                "store_entry_path": archive.name,
                "type": "bfabric_copy_resource",
            }
        ]
    }
    outputs_yml_path.write_text(yaml.safe_dump(document, sort_keys=False))
    return archive
