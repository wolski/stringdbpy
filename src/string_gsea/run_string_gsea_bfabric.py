#%% md

import shlex
import os
from pathlib import Path

import yaml
from loguru import logger

from string_gsea.gsea_utilities import get_rank_files, find_zip_files
from string_gsea.get_species import get_species_taxon
import subprocess
import tempfile

from string_gsea.postprocess import result_to_xlsx
from string_gsea.string_gsea_builder import StringGSEABuilder
from string_gsea.string_gsea_results import StringGSEAResults


def extract_workunit_id_from_file(file_path: Path) -> str | None:
    try:
        with open(file_path, 'r') as file:
            data = yaml.safe_load(file)
        # Navigate to the 'registration' block and extract 'workunit_id'
        return data.get('registration', {}).get('workunit_id')
    except Exception as e:
        logger.error(f"Error reading or parsing the YAML file: {e}")
        return None


def outputs_yml(search_zip : Path, base_path: Path, outputs_yml = "outputs.yml"):
    output1 = {
        'local_path': str(search_zip.resolve()),
        'store_entry_path': search_zip.name,
        'type': 'bfabric_copy_resource'
    }
    outputs_list = [output1]
    # Step 3: Create the main dictionary and assign the outputs list to the key 'outputs'
    data = {
        'outputs': outputs_list
    }
    outputs_yml = base_path / outputs_yml
    with open(outputs_yml, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)
    logger.info(f"YAML file {outputs_yml} has been generated.")


def _save_link(link:str, name:str, workunit_id: str) -> dict:
    # Define your arguments
    cmd = ["bfabric_save_link_to_workunit.py", str(workunit_id), link, name]
    logger.info(shlex.join(cmd))
    # Run the command and capture the output
    result = subprocess.run(cmd, capture_output=True, text=True)
    # Check if the command executed successfully
    if result.returncode == 0:
        return {result.returncode : result.stdout}
    else:
        return {result.returncode : result.stderr}


def save_link(links_dict: dict, workunit_id: str) -> dict:
    save_status = {}
    for outer_key, inner_dict in links_dict.items():
        for inner_key, link_url in inner_dict.items():
            link_name = f"Result string-db GSEA for [ {outer_key}/{inner_key} ]"
            save_status[(outer_key, inner_key)] = _save_link(link_url, link_name, workunit_id)
    return save_status


def register_result(workunit_id):
    # home/bfabric/slurmworker/bin/fgcz_app_runner 0.0.17 outputs register workunit_id
    # Define the command and its arguments as a list
    cmd = [
        "/home/bfabric/slurmworker/bin/fgcz_app_runner",
        "0.0.17",
        "outputs",
        "register",
        "outputs.yml",
        workunit_id
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result



def run_string_gsea_bfabric(
    zip_path: Path,
    workunit_id: str,
    fdr: float = 0.25,
    base_dir: Path = Path(".")
) -> None:
    # 1) Prepare config & workspace
    config = {
        "api_key": "b36F8oaRJwFZ",
        "fdr": fdr,
        "caller_identity": "www.fgcz.ch",
        "ge_enrichment_rank_direction": -1,
    }
    base_dir.mkdir(exist_ok=True)
    species = get_species_taxon(zip_path)
    dataframes = get_rank_files(zip_path)

    # 2) Build, write inputs, submit & poll
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

    # 3) Serialize session + results, write out files
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

    # 4) Post‑processing
    result_to_xlsx(tsv_dir, workunit_id)
    path = StringGSEAResults.zip_folder(results.get_res_path())
    logger.info(f"Zipped results to {path}")

    outputs_yml(path, base_dir)
    logger.info("Outputs YAML written.")

    # 5) Save link if desired
    if workunit_id:
        status = save_link(results.get_links(), workunit_id)
        logger.info(f"Saved link status: {status}")


if __name__ == "__main__":
    test_data = Path(__file__).parent.parent.parent / "tests/data/dummy_d"
    workunit_id = extract_workunit_id_from_file(test_data / "params.yml") or "876543"
    zip_path = test_data.parent / "2848501.zip"
    tempdir = Path(tempfile.mkdtemp())
    logger.info(f"Running on {zip_path} into {tempdir}")
    run_string_gsea_bfabric(zip_path, workunit_id, fdr=0.25, base_dir=tempdir)


