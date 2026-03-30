"""Centralized test fixtures for STRING-GSEA test suite."""

from pathlib import Path

import pytest

# Base paths
TEST_DIR = Path(__file__).parent
DATA_DIR = TEST_DIR / "data"
DATASETS_DIR = DATA_DIR / "datasets"
FIXTURES_DIR = DATA_DIR / "fixtures"
OUTPUTS_DIR = DATA_DIR / "outputs"


# ============================================================
# Dataset Fixtures (for Snakemake-compatible datasets)
# ============================================================


@pytest.fixture
def datasets_dir() -> Path:
    """Path to datasets directory containing Snakemake-iterable datasets."""
    return DATASETS_DIR


# ============================================================
# Input File Fixtures (backwards compatible)
# ============================================================


@pytest.fixture
def mouse_xlsx_zip() -> Path:
    """Mouse XLSX differential expression ZIP file."""
    return DATASETS_DIR / "mouse_xlsx" / "input.zip"


@pytest.fixture
def yeast_rnk_zip() -> Path:
    """Yeast RNK ZIP file."""
    return DATASETS_DIR / "yeast_rnk" / "input.zip"


@pytest.fixture
def human_rnk_2848501_zip() -> Path:
    """Human RNK ZIP file (2848501)."""
    return DATASETS_DIR / "human_rnk_2848501" / "input.zip"


# ============================================================
# Special Fixtures (test-specific, not for Snakemake)
# ============================================================


@pytest.fixture
def fasta_test_zip() -> Path:
    """FASTA test ZIP for species detection (no matching pattern)."""
    return FIXTURES_DIR / "no_matching_zip" / "fasta_test.zip"


@pytest.fixture
def no_matching_dir() -> Path:
    """Directory with non-matching ZIP files."""
    return FIXTURES_DIR / "no_matching_zip"


@pytest.fixture
def dummy_session_rnk_zip() -> Path:
    """Dummy session RNK ZIP file."""
    return FIXTURES_DIR / "dummy_session" / "DEA_mouse_fasta_rnk.zip"


@pytest.fixture
def dummy_session_yml() -> Path:
    """Dummy session YAML file."""
    return FIXTURES_DIR / "dummy_session" / "session.yml"


@pytest.fixture
def multi_contrast_tsv_dir() -> Path:
    """Directory with multiple contrast TSV files."""
    return OUTPUTS_DIR / "human_rnk_2848501" / "WU_2848501_GSEA" / "from_rnk"


@pytest.fixture
def single_contrast_tsv() -> Path:
    """Single-contrast STRING-DB GSEA TSV file for model tests."""
    return (
        OUTPUTS_DIR
        / "human_rnk_2848501"
        / "WU_2848501_GSEA"
        / "from_rnk"
        / "Bait_NCP_pUbT12_results.tsv"
    )


@pytest.fixture
def network_test_xlsx() -> Path:
    """XLSX file for network tests."""
    return (
        FIXTURES_DIR
        / "dummy_output"
        / "WU_abcd_GSEA"
        / "from_rnk"
        / "WUabcd_string_gsea_results_long.xlsx"
    )
