"""Integration tests for the full STRING-GSEA workflow (hits STRING API)."""

import shutil
import subprocess
from pathlib import Path
from typing import cast

import pytest

type DatasetParameters = tuple[str, str, str]

DATASETS_DIR = Path(__file__).parent / "data" / "datasets"
OUTPUTS_DIR = Path(__file__).parent / "data" / "outputs"

# (dataset_name, workunit_id, which)  — which="none" means RNK input
INTEGRATION_DATASETS: list[DatasetParameters] = [
    ("mouse_xlsx", "mouse_fasta", "pep_2_no_imputed"),
    ("human_rnk_2848501", "2848501", "none"),
    ("yeast_rnk", "yeast_rnk", "none"),
]

SMOKE_DATASETS: list[DatasetParameters] = [
    ("yeast_rnk", "yeast_rnk", "none"),
]


def _dataset_id(value: object) -> str:
    if not isinstance(value, tuple) or not value or not isinstance(value[0], str):
        raise TypeError("Dataset parameter must be a non-empty tuple beginning with a name")
    return value[0]


def _dataset_parameters(value: object) -> DatasetParameters:
    if not isinstance(value, tuple):
        raise TypeError("Dataset parameter must contain three strings")
    items = cast(tuple[object, ...], value)
    if len(items) != 3:
        raise TypeError("Dataset parameter must contain three strings")
    first, second, third = items
    if not isinstance(first, str) or not isinstance(second, str) or not isinstance(third, str):
        raise TypeError("Dataset parameter must contain three strings")
    return first, second, third


@pytest.fixture(params=INTEGRATION_DATASETS, ids=_dataset_id)
def dataset_params(request: pytest.FixtureRequest) -> DatasetParameters:
    return _dataset_parameters(request.param)


@pytest.fixture(params=SMOKE_DATASETS, ids=_dataset_id)
def smoke_params(request: pytest.FixtureRequest) -> DatasetParameters:
    return _dataset_parameters(request.param)


@pytest.mark.smoke
def test_workflow_smoke(smoke_params: DatasetParameters) -> None:
    """Quick single-contrast RNK workflow test (yeast, ~4 MB)."""
    _run_workflow(*smoke_params)


def _run_workflow(ds_name: str, workunit_id: str, which: str) -> None:
    """Run string_gsea_workflow end-to-end and verify outputs."""
    zip_path = str(DATASETS_DIR / ds_name / "input.zip")
    out_dir = OUTPUTS_DIR / ds_name

    # Clean previous run
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "string_gsea_workflow",
        zip_path,
        workunit_id,
        str(out_dir),
        "--which",
        which,
        "--cores",
        "1",
    ]
    result = subprocess.run(cmd)
    assert result.returncode == 0, f"Workflow failed for {ds_name}"

    # GSEA outputs
    wu_dir = out_dir / f"WU_{workunit_id}_GSEA"
    assert wu_dir.exists(), f"Expected output directory {wu_dir} not found"
    assert (wu_dir / f"WU{workunit_id}_gsea_result.json").exists()
    assert (wu_dir / "gsea_session.yml").exists()

    # Render outputs
    assert (out_dir / "render_done.txt").exists(), "Render step did not complete"

    # Package outputs
    assert (out_dir / "outputs.yml").exists(), "Package step did not produce outputs.yml"
    assert (out_dir / f"WU_{workunit_id}_GSEA.zip").exists(), "Package step did not produce zip"


@pytest.mark.integration
def test_workflow_end_to_end(dataset_params: DatasetParameters) -> None:
    """Run string_gsea_workflow end-to-end: GSEA -> render -> package."""
    _run_workflow(*dataset_params)
