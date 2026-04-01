from pathlib import Path
from typing import Literal

from cyclopts import App
from loguru import logger

from string_gsea.gsea_config import get_configuration
from string_gsea.gsea_result_processor import write_gsea_xlsx
from string_gsea.gsea_utilities import get_rank_files, write_rank_files
from string_gsea.models.gsea_models import RunMetadata, parse_gsea_results
from string_gsea.ranks_from_dea_xlsx import DiffXLSX
from string_gsea.species_detection import get_species_taxon
from string_gsea.string_gsea_builder import StringGSEABuilder
from string_gsea.string_gsea_results import StringGSEAResults

app = App()


@app.default()
def string_gsea_run(
    zip_path: str,
    workunit_id: str,
    out_dir: str = ".",
    from_rnk: bool = False,
    which: Literal["pep_1", "pep_1_no_imputed", "pep_2", "pep_2_no_imputed"] = "pep_2_no_imputed",
    create_zip: bool = False,
):
    """
    Run STRING GSEA analysis on the provided zip file.

    Args:
        zip_path (str): Path to the input zip file containing rank data.
        workunit_id (str): Identifier for this analysis run.
        out_dir (str, optional): Base directory for output files. Defaults to ".".
    """
    zip_path = Path(zip_path)
    base_dir = Path(out_dir)

    # 1) Configuration
    config = get_configuration()
    base_dir.mkdir(exist_ok=True)
    if not zip_path.exists():
        raise FileNotFoundError(f"Zip file not found: {zip_path}")

    # 2) Read rank data
    if from_rnk:
        rank_lists = get_rank_files(zip_path)
    else:
        df_xlsx = DiffXLSX(zip_path)
        rank_lists = df_xlsx.rank_dict(which=which)

    # 3) Detect species
    species = get_species_taxon(zip_path, rank_lists, api_base_url=config.api_base_url)

    # 4) Submit & poll
    builder = StringGSEABuilder(
        rank_lists=rank_lists,
        config=config,
        workunit_id=workunit_id,
        species=species,
        base_path=base_dir,
    )
    builder.submit().poll()

    # 5) Output directory
    res_path = base_dir / f"WU_{workunit_id}_GSEA"
    res_path.mkdir(parents=True, exist_ok=True)

    session_yaml = builder.save_session(res_path / "gsea_session.yml")
    logger.info(f"Session YAML written to {session_yaml}")

    # 6) Download raw results into memory
    raw_results = StringGSEAResults(builder.session)
    raw_results.download()

    # 7) Parse into typed model
    metadata = RunMetadata.from_config(config, workunit_id=workunit_id, species=species)
    gsea_result = parse_gsea_results(rank_lists, raw_results.tsv_content, metadata=metadata)

    # 8) Write outputs
    write_rank_files(rank_lists, res_path)
    gsea_result.to_json(res_path / f"WU{workunit_id}_gsea_result.json")
    write_gsea_xlsx(gsea_result, workunit_id, res_path)

    # 9) Persist raw files to disk
    raw_results.write_tsv(res_path)
    raw_results.write_graphs(res_path)
    raw_results.write_links(res_path)

    if create_zip:
        path = StringGSEAResults.zip_folder(res_path)
        logger.info(f"Zipped results to {path}")


if __name__ == "__main__":
    app()
