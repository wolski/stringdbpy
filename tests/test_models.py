import pytest

from string_gsea.models.gsea_models import (
    CategoryGSEA,
    GenePool,
    GSEAResult,
    TermGSEA,
    parse_gsea_tsv,
    parse_gsea_tsv_dir,
    parse_rank_file,
)


@pytest.fixture
def parsed_categories(single_contrast_tsv):
    return parse_gsea_tsv(single_contrast_tsv)


def test_parse_returns_all_categories(parsed_categories):
    assert len(parsed_categories) > 10
    assert "KEGG" in parsed_categories
    assert "GO Process" in parsed_categories


def test_parse_category_filter(single_contrast_tsv):
    result = parse_gsea_tsv(single_contrast_tsv, categories={"KEGG"})
    assert list(result.keys()) == ["KEGG"]


def test_gene_pool_shared_across_categories(parsed_categories):
    """All categories for one contrast share the same gene pool object."""
    pools = [cat.gene_pool for cat in parsed_categories.values()]
    assert all(p is pools[0] for p in pools)


def test_gene_hit_values(parsed_categories):
    """Check a known gene from the first row (GO:0006364)."""
    cat = parsed_categories["GO Process"]
    term = next(t for t in cat.terms if t.term_id == "GO:0006364")
    first_pid = term.gene_ids[0]
    hit = cat.gene_pool[first_pid]
    assert hit.protein_id == "9606.ENSP00000007390"
    assert hit.label == "TSR3"
    assert hit.input_label == "Q9UJK0"
    assert pytest.approx(hit.input_value, rel=1e-3) == 1.8687
    assert hit.rank == 2241


def test_term_gsea_values(parsed_categories):
    cat = parsed_categories["GO Process"]
    term = next(t for t in cat.terms if t.term_id == "GO:0006364")
    assert pytest.approx(term.enrichment_score, rel=1e-4) == 2.09208
    assert term.direction == "bottom"
    assert pytest.approx(term.fdr, rel=1e-2) == 6.99e-20
    assert term.genes_mapped == 168
    assert term.genes_in_set == 220
    assert len(term.gene_ids) == 168


def test_contrast_from_filename(parsed_categories):
    cat = next(iter(parsed_categories.values()))
    assert cat.contrast == "Bait_NCP_pUbT12_results.tsv"


def test_contrast_override(single_contrast_tsv):
    result = parse_gsea_tsv(single_contrast_tsv, contrast="custom_name", categories={"KEGG"})
    assert result["KEGG"].contrast == "custom_name"


def test_frozen_immutability(parsed_categories):
    cat = next(iter(parsed_categories.values()))
    hit = next(iter(cat.gene_pool.values()))
    with pytest.raises(AttributeError):
        hit.label = "changed"
    term = cat.terms[0]
    with pytest.raises(AttributeError):
        term.fdr = 0.99


def test_invalid_pool_reference():
    pool = GenePool(entries={})
    term = TermGSEA(
        term_id="GO:0000001",
        category="test",
        description="test",
        enrichment_score=1.0,
        direction="top",
        fdr=0.01,
        method="ks",
        genes_mapped=1,
        genes_in_set=10,
        gene_ids=("nonexistent_protein",),
    )
    with pytest.raises(ValueError, match="unknown protein_ids"):
        CategoryGSEA(
            category="test",
            contrast="test",
            gene_pool=pool,
            terms=(term,),
        )


# ---------------------------------------------------------------------------
# Container + multi-contrast tests
# ---------------------------------------------------------------------------


@pytest.fixture
def gsea_result(multi_contrast_tsv_dir) -> GSEAResult:
    return parse_gsea_tsv_dir(multi_contrast_tsv_dir)


def test_parse_dir_returns_both_contrasts(gsea_result):
    assert len(gsea_result.contrast_names) == 2
    assert any("pUbT12_" in c for c in gsea_result.contrast_names)
    assert any("pUbT12T14_" in c for c in gsea_result.contrast_names)


def test_gsea_result_category_names(gsea_result):
    cats = gsea_result.category_names
    assert "KEGG" in cats
    assert "GO Process" in cats
    assert len(cats) > 10


def test_get_multi_category(gsea_result):
    contrast = gsea_result.contrast_names[0]
    mc = gsea_result.get_multi_category(contrast)
    assert mc.contrast == contrast
    assert "KEGG" in mc.categories


def test_get_multi_contrast(gsea_result):
    mc = gsea_result.get_multi_contrast("KEGG")
    assert mc.category == "KEGG"
    assert len(mc.contrasts) == 2


def test_get_category(gsea_result):
    contrast = gsea_result.contrast_names[0]
    cat = gsea_result.get_category(contrast, "KEGG")
    assert isinstance(cat, CategoryGSEA)
    assert cat.category == "KEGG"
    assert cat.contrast == contrast


# ---------------------------------------------------------------------------
# Computed scores tests
# ---------------------------------------------------------------------------


def test_gene_pool_n_genes(parsed_categories):
    """Shared gene pool has n_genes property."""
    cat = parsed_categories["GO Process"]
    assert cat.gene_pool.n_genes > 1000
    assert cat.gene_pool.n_genes == len(cat.gene_pool)


def test_gene_ratio(parsed_categories):
    cat = parsed_categories["GO Process"]
    term = next(t for t in cat.terms if t.term_id == "GO:0006364")
    assert pytest.approx(term.gene_ratio, rel=1e-2) == 168 / 220


def test_mean_input_value(parsed_categories):
    cat = parsed_categories["GO Process"]
    term = next(t for t in cat.terms if t.term_id == "GO:0006364")
    miv = term.mean_input_value(cat.gene_pool)
    assert isinstance(miv, float)
    assert miv != 0.0


def test_rank_nes(multi_contrast_tsv_dir):
    """rank_nes needs total ranked genes (from RankList), not per-category pool size."""
    rl = parse_rank_file(multi_contrast_tsv_dir / "Bait_NCP_pUbT12.rnk")
    cats = parse_gsea_tsv(
        multi_contrast_tsv_dir / "Bait_NCP_pUbT12_results.tsv",
        categories={"GO Process"},
    )
    cat = cats["GO Process"]
    term = next(t for t in cat.terms if t.term_id == "GO:0006364")
    nes = term.rank_nes(cat.gene_pool, rl.n_genes)
    # Also verify it works with gene_pool.n_genes (mapped count)
    nes_pool = term.rank_nes(cat.gene_pool, cat.gene_pool.n_genes)
    assert isinstance(nes_pool, float)
    assert -1.0 <= nes <= 1.0
    # direction is "bottom" → genes are at the bottom of the list → negative NES
    assert nes < 0


# ---------------------------------------------------------------------------
# RankList + mapping efficiency tests
# ---------------------------------------------------------------------------


def test_parse_rank_file(multi_contrast_tsv_dir):
    rnk_path = multi_contrast_tsv_dir / "Bait_NCP_pUbT12.rnk"
    rl = parse_rank_file(rnk_path)
    assert rl.contrast == "Bait_NCP_pUbT12"
    assert rl.n_genes > 3000
    assert "A0AV96" in rl.entries
    assert pytest.approx(rl.entries["A0AV96"], rel=1e-3) == 2.2147


def test_gsea_result_has_rank_lists(gsea_result):
    assert len(gsea_result.rank_lists) == 2
    for contrast in gsea_result.contrast_names:
        assert contrast in gsea_result.rank_lists


def test_mapping_efficiency(gsea_result):
    for contrast in gsea_result.contrast_names:
        eff = gsea_result.mapping_efficiency(contrast)
        assert 0.5 < eff < 1.0
