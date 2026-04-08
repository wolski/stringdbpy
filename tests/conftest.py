"""Centralized test fixtures for STRING-GSEA test suite."""

from pathlib import Path

import pytest

from string_gsea.gsea_config import GSEAConfig

# Base paths
TEST_DIR = Path(__file__).parent
DATA_DIR = TEST_DIR / "data"
DATASETS_DIR = DATA_DIR / "datasets"
FIXTURES_DIR = DATA_DIR / "fixtures"
OUTPUTS_DIR = DATA_DIR / "outputs"


# ============================================================
# Dataset Fixtures (for Snakemake-compatible datasets)
# ============================================================


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


# ============================================================
# Special Fixtures (test-specific, not for Snakemake)
# ============================================================


@pytest.fixture
def fasta_test_zip() -> Path:
    """FASTA test ZIP for species detection (no matching pattern)."""
    return FIXTURES_DIR / "no_matching_zip" / "fasta_test.zip"


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
    d = OUTPUTS_DIR / "human_rnk_2848501" / "WU_2848501_GSEA" / "from_rnk"
    if not d.exists():
        pytest.skip("integration output not available (run 'make test-integration' first)")
    return d


@pytest.fixture
def single_contrast_tsv() -> Path:
    """Single-contrast STRING-DB GSEA TSV file for model tests."""
    f = OUTPUTS_DIR / "human_rnk_2848501" / "WU_2848501_GSEA" / "from_rnk" / "Bait_NCP_pUbT12_results.tsv"
    if not f.exists():
        pytest.skip("integration output not available (run 'make test-integration' first)")
    return f


@pytest.fixture
def gsea_config() -> GSEAConfig:
    """Minimal GSEAConfig for unit tests (no real API key)."""
    return GSEAConfig(
        api_key="test_key",
        fdr=0.25,
        ge_enrichment_rank_direction=1,
        caller_identity="test@example.com",
    )
