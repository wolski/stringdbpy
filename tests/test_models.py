import polars as pl
import pytest

from string_gsea.models.gsea_models import (
    CategoryGSEA,
    GenePool,
    GSEAResult,
    RankList,
    RankListCollection,
    RunMetadata,
    TermGSEA,
    parse_gsea_results,
    parse_gsea_tsv,
    parse_gsea_tsv_dir,
    parse_gsea_tsv_from_string,
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


# ---------------------------------------------------------------------------
# to_polars_long tests
# ---------------------------------------------------------------------------


def test_to_polars_long_columns(gsea_result):
    df = gsea_result.to_polars_long()
    expected_cols = {
        "contrast",
        "category",
        "termID",
        "termDescription",
        "enrichmentScore",
        "direction",
        "falseDiscoveryRate",
        "method",
        "genesMapped",
        "genesInSet",
        "proteinIDs",
        "proteinLabels",
        "proteinInputLabels",
        "proteinInputValues",
        "proteinRanks",
        "directionNR",
        "num_contrasts",
    }
    assert set(df.columns) == expected_cols


def test_to_polars_long_row_count(gsea_result):
    df = gsea_result.to_polars_long()
    # Row count should equal total number of terms across all contrasts/categories
    total_terms = sum(len(cat.terms) for mc in gsea_result.data.values() for cat in mc.categories.values())
    assert len(df) == total_terms


def test_to_polars_long_values(gsea_result):
    """Spot-check a known term's values in the long DataFrame."""
    df = gsea_result.to_polars_long()
    # Filter for a known term from the first contrast
    # Use the pUbT12 contrast (not pUbT12T14) where we know exact values
    contrast = next(c for c in gsea_result.contrast_names if "pUbT12_" in c)
    row = df.filter((pl.col("contrast") == contrast) & (pl.col("termID") == "GO:0006364"))
    assert len(row) == 1
    assert row["category"][0] == "GO Process"
    assert row["direction"][0] == "bottom"
    assert row["directionNR"][0] == -1
    assert row["genesMapped"][0] == 168
    assert row["genesInSet"][0] == 220
    # Protein columns should be comma-separated strings
    assert "," in row["proteinIDs"][0]


# ---------------------------------------------------------------------------
# JSON serialization tests
# ---------------------------------------------------------------------------


def test_json_round_trip(gsea_result, tmp_path):
    """Serialize to JSON and deserialize — result should be identical."""
    json_path = tmp_path / "gsea_result.json"
    gsea_result.to_json(json_path)
    restored = GSEAResult.from_json(json_path)

    assert restored.contrast_names == gsea_result.contrast_names
    assert restored.category_names == gsea_result.category_names
    assert len(restored.rank_lists) == len(gsea_result.rank_lists)

    # Spot-check a term
    contrast = gsea_result.contrast_names[0]
    cat_orig = gsea_result.get_category(contrast, "KEGG")
    cat_restored = restored.get_category(contrast, "KEGG")
    assert len(cat_restored.terms) == len(cat_orig.terms)
    assert cat_restored.terms[0].term_id == cat_orig.terms[0].term_id
    assert cat_restored.terms[0].fdr == cat_orig.terms[0].fdr
    assert cat_restored.terms[0].gene_ids == cat_orig.terms[0].gene_ids


def test_json_shared_gene_pool(gsea_result, tmp_path):
    """After round-trip, gene pool is still shared across categories."""
    json_path = tmp_path / "gsea_result.json"
    gsea_result.to_json(json_path)
    restored = GSEAResult.from_json(json_path)

    contrast = restored.contrast_names[0]
    mc = restored.get_multi_category(contrast)
    pools = [cat.gene_pool for cat in mc.categories.values()]
    assert all(p is pools[0] for p in pools)


def test_json_gene_hit_values(gsea_result, tmp_path):
    """Gene hit values survive round-trip."""
    json_path = tmp_path / "gsea_result.json"
    gsea_result.to_json(json_path)
    restored = GSEAResult.from_json(json_path)

    # Use the pUbT12 contrast
    contrast = next(c for c in restored.contrast_names if "pUbT12_" in c)
    cat = restored.get_category(contrast, "GO Process")
    term = next(t for t in cat.terms if t.term_id == "GO:0006364")
    hit = cat.gene_pool[term.gene_ids[0]]
    assert hit.protein_id == "9606.ENSP00000007390"
    assert hit.label == "TSR3"
    assert pytest.approx(hit.input_value, rel=1e-3) == 1.8687


def test_json_rank_list_round_trip(gsea_result, tmp_path):
    """Rank lists survive round-trip."""
    json_path = tmp_path / "gsea_result.json"
    gsea_result.to_json(json_path)
    restored = GSEAResult.from_json(json_path)

    for contrast in gsea_result.contrast_names:
        orig_rl = gsea_result.rank_lists[contrast]
        rest_rl = restored.rank_lists[contrast]
        assert rest_rl.contrast == orig_rl.contrast
        assert rest_rl.n_genes == orig_rl.n_genes


# ---------------------------------------------------------------------------
# RankList — new features
# ---------------------------------------------------------------------------


def test_rank_list_to_rnk_string():
    rl = RankList(contrast="c", entries={"GENE1": 1.5, "GENE2": -0.3})
    text = rl.to_rnk_string()
    lines = text.strip().split("\n")
    assert len(lines) == 2
    assert lines[0] == "GENE1\t1.5"
    assert lines[1] == "GENE2\t-0.3"


def test_rank_list_sample_identifiers():
    entries = {f"G{i}": float(i) for i in range(20)}
    rl = RankList(contrast="c", entries=entries)
    sampled = rl.sample_identifiers(5)
    assert len(sampled) == 5
    assert all(s in entries for s in sampled)


def test_rank_list_sample_identifiers_fewer_than_nr():
    rl = RankList(contrast="c", entries={"G1": 1.0, "G2": 2.0})
    sampled = rl.sample_identifiers(10)
    assert len(sampled) == 2


def test_rank_list_from_polars():
    df = pl.DataFrame({"id": ["GENE1", "GENE2"], "statistic": [1.5, -0.3]})
    rl = RankList.from_polars(df, contrast="test")
    assert rl.contrast == "test"
    assert rl.entries == {"GENE1": 1.5, "GENE2": -0.3}


def test_rank_list_json_round_trip():
    rl = RankList(contrast="c", entries={"G1": 1.0})
    d = rl.to_dict()
    restored = RankList.from_dict(d)
    assert restored == rl


# ---------------------------------------------------------------------------
# RankListCollection
# ---------------------------------------------------------------------------


def test_rank_list_collection_basic():
    rl1 = RankList(contrast="A", entries={"G1": 1.0})
    rl2 = RankList(contrast="B", entries={"G2": 2.0})
    coll = RankListCollection(analysis="pep_1", rank_lists=[rl1, rl2])
    assert coll.analysis == "pep_1"
    assert len(coll) == 2
    assert coll["A"] is rl1
    assert "B" in coll


def test_rank_list_collection_iteration():
    rl1 = RankList(contrast="A", entries={"G1": 1.0})
    rl2 = RankList(contrast="B", entries={"G2": 2.0})
    coll = RankListCollection(analysis="pep_1", rank_lists=[rl1, rl2])
    items = list(coll)
    assert len(items) == 2
    assert rl1 in items
    assert rl2 in items


def test_rank_list_collection_contrasts():
    rl1 = RankList(contrast="A", entries={"G1": 1.0})
    rl2 = RankList(contrast="B", entries={"G2": 2.0})
    coll = RankListCollection(analysis="pep_2", rank_lists=[rl1, rl2])
    assert coll.contrasts == ["A", "B"]


def test_rank_list_collection_first():
    rl1 = RankList(contrast="A", entries={"G1": 1.0})
    coll = RankListCollection(analysis="pep_1", rank_lists=[rl1])
    assert coll.first() is rl1


def test_rank_list_collection_add():
    coll = RankListCollection(analysis="pep_1", rank_lists=[])
    assert len(coll) == 0
    coll.add(RankList(contrast="A", entries={"G1": 1.0}))
    assert len(coll) == 1


# ---------------------------------------------------------------------------
# parse_gsea_tsv_from_string tests
# ---------------------------------------------------------------------------


def test_parse_gsea_tsv_from_string(single_contrast_tsv):
    """parse_gsea_tsv_from_string should produce the same result as parse_gsea_tsv."""
    content = single_contrast_tsv.read_text()
    contrast = single_contrast_tsv.name

    from_file = parse_gsea_tsv(single_contrast_tsv, contrast=contrast)
    from_string = parse_gsea_tsv_from_string(content, contrast=contrast)

    assert set(from_file.keys()) == set(from_string.keys())
    for cat_name in from_file:
        assert len(from_file[cat_name].terms) == len(from_string[cat_name].terms)


def test_parse_gsea_tsv_from_string_category_filter(single_contrast_tsv):
    content = single_contrast_tsv.read_text()
    result = parse_gsea_tsv_from_string(content, contrast="test", categories={"KEGG"})
    assert list(result.keys()) == ["KEGG"]


# ---------------------------------------------------------------------------
# parse_gsea_results tests
# ---------------------------------------------------------------------------


def test_parse_gsea_results(single_contrast_tsv):
    """parse_gsea_results should build GSEAResult from in-memory TSV content."""
    content = single_contrast_tsv.read_text()

    rank_lists = RankListCollection(
        analysis="from_rnk",
        rank_lists=[RankList(contrast="Bait_NCP_pUbT12", entries={"A0AV96": 2.2147})],
    )
    tsv_content = {("from_rnk", "Bait_NCP_pUbT12"): content}

    gsea_result = parse_gsea_results(rank_lists, tsv_content)

    assert len(gsea_result.contrast_names) == 1
    assert "Bait_NCP_pUbT12_results.tsv" in gsea_result.contrast_names
    assert "KEGG" in gsea_result.category_names
    assert "Bait_NCP_pUbT12_results.tsv" in gsea_result.rank_lists
    assert gsea_result.metadata is None


def test_run_metadata_round_trip():
    """RunMetadata should survive dict serialization round-trip."""
    meta = RunMetadata(
        workunit_id="WU123",
        species=9606,
        fdr=0.25,
        ge_enrichment_rank_direction=1,
        caller_identity="test@example.com",
        api_base_url="https://version-12-0.string-db.org/api",
    )
    d = meta.to_dict()
    restored = RunMetadata.from_dict(d)
    assert restored == meta


def test_run_metadata_from_config():
    """RunMetadata.from_config should extract non-sensitive fields from GSEAConfig."""
    from string_gsea.gsea_config import GSEAConfig

    config = GSEAConfig(
        api_key="secret-key-123",
        fdr=0.05,
        ge_enrichment_rank_direction=1,
        caller_identity="test@example.com",
        api_base_url="https://version-12-0.string-db.org/api",
    )
    meta = RunMetadata.from_config(config, workunit_id="WU456", species=10090)
    assert meta.workunit_id == "WU456"
    assert meta.species == 10090
    assert meta.fdr == 0.05
    assert meta.caller_identity == "test@example.com"
    # api_key must NOT leak into RunMetadata
    assert not hasattr(meta, "api_key")


def test_gsea_result_with_metadata_round_trip(single_contrast_tsv):
    """GSEAResult with metadata should survive JSON round-trip."""
    content = single_contrast_tsv.read_text()
    rank_lists = RankListCollection(
        analysis="from_rnk",
        rank_lists=[RankList(contrast="Bait_NCP_pUbT12", entries={"A0AV96": 2.2147})],
    )
    tsv_content = {("from_rnk", "Bait_NCP_pUbT12"): content}
    meta = RunMetadata(
        workunit_id="WU789",
        species=9606,
        fdr=0.25,
        ge_enrichment_rank_direction=1,
        caller_identity="test@example.com",
        api_base_url="https://version-12-0.string-db.org/api",
    )

    gsea_result = parse_gsea_results(rank_lists, tsv_content, metadata=meta)
    assert gsea_result.metadata == meta

    # JSON round-trip
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        json_path = Path(tmp) / "result.json"
        gsea_result.to_json(json_path)
        restored = GSEAResult.from_json(json_path)

    assert restored.metadata is not None
    assert restored.metadata.workunit_id == "WU789"
    assert restored.metadata.species == 9606
    assert restored.metadata.fdr == 0.25
