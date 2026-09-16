# pyright: reportUnknownMemberType=false, reportMissingTypeStubs=false
"""Boundary tests for the AnnData DEA-artifact reader.

The fixtures in `test_gsea_anndata_input.py` cover the artifacts prolfquapp
actually writes; these cover what the reader does with one that is malformed,
which is where a silent wrong answer would hurt most.
"""

from collections.abc import Mapping
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from numpy.typing import NDArray

from string_gsea.gsea.dea_artifact import ContrastRoles, read_dea_artifact

LM_ROLES: dict[str, object] = {
    "subject_id": "protein_Id",
    "contrast_col": "contrast",
    "effect_col": "diff",
    "score_col": "statistic",
    "pvalue_col": "p.value",
    "fdr_col": "FDR",
    "model_name_col": "modelName",
    "significance_directional": False,
}


def _write(
    path: Path,
    metadata: dict[str, object],
    *,
    varm: Mapping[str, NDArray[np.float64]] | None = None,
) -> Path:
    """Write a minimal two-feature artifact carrying the given metadata."""
    import pandas as pd

    adata = ad.AnnData(
        X=np.zeros((2, 2), dtype=float),
        obs=pd.DataFrame(index=["S1", "S2"]),
        var=pd.DataFrame({"protein_Id": ["P1", "P2"]}, index=["P1", "P2"]),
    )
    if varm is not None:
        for key, values in varm.items():
            adata.varm[key] = values
    adata.uns["prolfquapp"] = metadata
    adata.write_h5ad(path)
    return path


def _metadata(**overrides: object) -> dict[str, object]:
    metadata: dict[str, object] = {
        "artifact_type": "dea_results",
        "identifier_key": "protein_Id",
        "contrast_configuration": dict(LM_ROLES),
        "varm_key_order": ["constrast_A"],
        "varm_columns": {"constrast_A": ["diff", "statistic", "p.value", "FDR"]},
        "varm_annotations": {"constrast_A": {"contrast": ["A", "A"], "modelName": ["lm", "lm"]}},
    }
    metadata.update(overrides)
    return metadata


def _values() -> dict[str, NDArray[np.float64]]:
    return {
        "constrast_A": np.array(
            [[1.0, 2.0, 0.01, 0.02], [-1.0, -2.0, 0.2, 0.3]],
            dtype=float,
        )
    }


def test_reads_a_minimal_artifact(tmp_path: Path) -> None:
    path = _write(tmp_path / "ok.h5ad", _metadata(), varm=_values())

    artifact = read_dea_artifact(path)

    assert artifact.identifier == "protein_Id"
    assert artifact.contrasts == ["A"]
    assert artifact.roles.tests_differences
    assert artifact.rows.height == 2


def test_wrong_artifact_type_is_refused(tmp_path: Path) -> None:
    path = _write(tmp_path / "lfq.h5ad", _metadata(artifact_type="lfqdata"), varm=_values())

    with pytest.raises(ValueError, match="not a prolfquapp dea_results artifact"):
        read_dea_artifact(path)


def test_missing_column_roles_are_refused(tmp_path: Path) -> None:
    metadata = _metadata()
    del metadata["contrast_configuration"]
    path = _write(tmp_path / "no-roles.h5ad", metadata, varm=_values())

    with pytest.raises(ValueError, match="no contrast_configuration"):
        read_dea_artifact(path)


def test_identifier_key_must_name_an_annotation_column(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "bad-id.h5ad",
        _metadata(identifier_key="does_not_exist"),
        varm=_values(),
    )

    with pytest.raises(ValueError, match="identifier_key names 'does_not_exist'"):
        read_dea_artifact(path)


def test_varm_column_names_must_match_the_matrix(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "short.h5ad",
        _metadata(varm_columns={"constrast_A": ["diff", "statistic"]}),
        varm=_values(),
    )

    with pytest.raises(ValueError, match="varm_columns names 2"):
        read_dea_artifact(path)


def test_an_artifact_without_contrasts_is_refused(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "empty.h5ad",
        _metadata(varm_key_order=["stats_raw_wide"]),
        varm=_values(),
    )

    with pytest.raises(ValueError, match="carries no contrast results"):
        read_dea_artifact(path)


def test_missing_varm_key_order_is_refused(tmp_path: Path) -> None:
    metadata = _metadata()
    del metadata["varm_key_order"]
    path = _write(tmp_path / "unordered.h5ad", metadata, varm=_values())

    with pytest.raises(ValueError, match="Missing varm_key_order"):
        read_dea_artifact(path)


def test_roles_without_a_pvalue_rank_on_the_effect() -> None:
    roles = ContrastRoles.from_metadata(
        {"contrast_configuration": {**LM_ROLES, "pvalue_col": "NA"}}
    )

    assert not roles.tests_differences
    assert roles.pvalue_col is None
    assert "diff" in str(roles.rank_score())


def test_roles_read_values_however_hdf5_boxed_them() -> None:
    roles = ContrastRoles.from_metadata(
        {
            "contrast_configuration": {
                **LM_ROLES,
                "subject_id": np.array(["protein_Id"]),
                "contrast_col": b"contrast",
                "significance_directional": np.array([True]),
            }
        }
    )

    assert roles.subject_id == "protein_Id"
    assert roles.contrast_col == "contrast"
    assert roles.directional


def test_a_role_that_is_not_text_is_refused() -> None:
    with pytest.raises(ValueError, match="Expected a string for subject_id"):
        ContrastRoles.from_metadata({"contrast_configuration": {**LM_ROLES, "subject_id": 42}})
