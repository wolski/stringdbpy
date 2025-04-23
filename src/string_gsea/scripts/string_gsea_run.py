from cyclopts import App
from pathlib import Path
from loguru import logger

from string_gsea.config import get_configuration
from string_gsea.gsea_utilities import get_rank_files
from string_gsea.get_species import get_species_taxon
from string_gsea.string_gsea_builder import StringGSEABuilder
from string_gsea.string_gsea_results import StringGSEAResults
from string_gsea.postprocess import result_to_xlsx

app = App()

@app.default()
def string_gsea_run(
    zip_path: str,
    workunit_id: str,
    out_dir: str = ".",
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
    dataframes = get_rank_files(zip_path)
    
    # 3) Build, write inputs, submit & poll
    builder = StringGSEABuilder(
        rank_dataframes=dataframes,
        config_dict=config,
        workunit_id=workunit_id,
        species=species,
        base_path=base_dir
    )
    builder.write_rank_files()
    # submit + poll under the hood
    results: StringGSEAResults = builder.get_result()
    logger.info(f"Jobs completed: {builder.session.res_job_id}")

    # 4) Serialize session + results, write out files
    #   - YAML session
    session_yaml = builder.save_session()
    logger.info(f"Session YAML written to {session_yaml}")

    #   - JSON results
    serialized_json = results.serialize_results()
    logger.info(f"Results JSON written to {serialized_json}")

    #   - links, TSVs, graphs
    links = results.write_links()
    tsv_dir = results.write_gsea_tsv()
    graph_dir = results.write_gsea_graphs()
    logger.info(f"Wrote links to {links}\nTSVs to {tsv_dir}\nGraphs to {graph_dir}")

    # 5) Post‑processing
    result_to_xlsx(tsv_dir, workunit_id)
    path = StringGSEAResults.zip_folder(results.get_res_path())
    logger.info(f"Zipped results to {path}")

if __name__ == '__main__':
    app()
