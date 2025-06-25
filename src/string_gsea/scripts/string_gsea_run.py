from cyclopts import App
from pathlib import Path
from loguru import logger

from string_gsea.gsea_config import get_configuration
from string_gsea.gsea_utilities import get_rank_files
from string_gsea.get_species import get_species_taxon
from string_gsea.string_gsea_builder import StringGSEABuilder
from string_gsea.string_gsea_results import StringGSEAResults
from string_gsea.gsea_result_processor import GSEAResultProcessor
from string_gsea.ranks_from_dea_xlsx import DiffXLSX

app = App()

@app.default()
def string_gsea_run(
    zip_path: str,
    workunit_id: str,
    out_dir: str = ".",
    from_rnk: bool = False,
    which: str = "pep_2_no_imputed"
):
    """
    Run STRING GSEA analysis on the provided zip file.

    Args:
        zip_path (str): Path to the input zip file containing rank data.
        workunit_id (str): Identifier for this analysis run.
        out_dir (str, optional): Base directory for output files. Defaults to ".".
    """
    # Convert paths to Path objects
    zip_path = Path(zip_path)
    base_dir = Path(out_dir)
    
    # 1) Get configuration from config file
    config = get_configuration()
    base_dir.mkdir(exist_ok=True)
    
    # 2) Get species and rank data
    species = get_species_taxon(zip_path)

    if from_rnk:
        dataframes = get_rank_files(zip_path)
    else:
        df_xlsx = DiffXLSX(zip_path)
        rank_files = df_xlsx.rank_dict(which=which)
        dataframes = rank_files
    
    # 3) Build, write inputs, submit & poll
    builder = StringGSEABuilder(
        rank_dataframes=dataframes,
        config=config,
        workunit_id=workunit_id,
        species=species,
        base_path=base_dir
    )
    xd = builder.get_res_path()
    print(xd)
    builder.write_rank_files()
    # submit + poll under the hood
    results: StringGSEAResults = builder.get_result()
    logger.info(f"Jobs completed: {builder.session.res_job_id}")

    # 4) Serialize session + results, write out files
    #   - YAML session
    session_yaml = builder.save_session()
    logger.info(f"Session YAML written to {session_yaml}")

    #   - JSON results
    # serialized_json = results.save_session()
    # logger.info(f"Results JSON written to {serialized_json}")

    #   - links, TSVs, graphs
    links = results.write_links()
    tsv_dir = results.write_gsea_tsv()
    graph_dir = results.write_gsea_graphs()
    logger.info(f"Wrote links to {links}\nTSVs to {tsv_dir}\nGraphs to {graph_dir}")

    # 5) Post‑processing
    GSEAResultProcessor.result_to_xlsx(tsv_dir, workunit_id)
    path = StringGSEAResults.zip_folder(results.get_res_path())
    logger.info(f"Zipped results to {path}")

def test_run():
    # Test code
    wd = Path(__file__).parent.parent.parent.parent
    zip_path = wd / "tests/data/DE_mouse_fasta_xlsx.zip"
    base_dir = Path(wd/"tests/data/dummy_res_2")
    workunit_id = 1234
    which = "pep_1"
    from_rnk = False

    logger.info("Running test with parameters:")
    logger.info(f"Zip path: {zip_path}")
    logger.info(f"Base directory: {base_dir}")
    logger.info(f"Workunit ID: {workunit_id}")
    logger.info(f"Analysis type: {which}")
    logger.info(f"From RNK: {from_rnk}")

    string_gsea_run(
        zip_path=zip_path,
        workunit_id=workunit_id,
        out_dir=base_dir,
        which=which,
        from_rnk=from_rnk
    )

if __name__ == '__main__':
    # test_run()
    # else:
    app()
