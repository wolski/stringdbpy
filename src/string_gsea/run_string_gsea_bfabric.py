#%% md

import shlex
import os
from pathlib import Path

import yaml
from loguru import logger

from string_gsea import postprocess
from string_gsea.gsea_utilities import get_rank_files, find_zip_files
from string_gsea.string_gsea import StringGSEA
from string_gsea.get_species import OxFieldsZip, get_species_taxon
import subprocess
import tempfile


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


def run_string_gsea_bfabric(zip_path : Path,
                            workunit_id: str,
                            fdr: float = 0.25,
                            base_dir: Path = Path(".")) -> None:
    api_key = "b36F8oaRJwFZ"
    base_dir.mkdir(exist_ok = True)
    species = get_species_taxon(zip_path)
    dataframes = get_rank_files(zip_path)
    gsea = StringGSEA(api_key, workunit_id, dataframes, species, fdr, base_dir)
    gsea.submit()
    logger.info(f"Job submitted successfully.{gsea.res_job_id}")
    gsea.pull_results()
    gsea.serialize_results()

    logger.info("got results")
    gsea.write_rank_files()
    result_dir = gsea.write_gsea_results()
    written_links = gsea.write_links()

    postprocess.result_to_xlsx(result_dir, workunit_id)
    logger.info(f"Results written to {result_dir}")

    path = gsea.zip_folder(result_dir)
    logger.info(f"Zipped results to {path}")
    outputs_yml(path, base_dir)
    logger.info("zipped results.")
    if workunit_id is not None:
        status = save_link(gsea.get_links(), workunit_id)
        logger.info(f"saved link {status}")


if __name__ == '__main__':
    test_data = Path(__file__).parent.parent.parent / "tests/data/"
    testing = True
    fdr: float = 0.25
    current_directory = os.getcwd()
    # Print the working directory
    print("Current working directory:", current_directory)
    workunit_id =  extract_workunit_id_from_file(test_data / "params.yml")
    if workunit_id is None:
        workunit_id = "876543"

    logger.info(f"Workunit ID: {workunit_id}")
    if testing:
        # replace the code above with pathlib and __FILE__
        zip_path = test_data / "DE_mouse_fasta_rnk.zip"
        logger.info(f"Zip path: {zip_path}")
    else:
        zip_path = find_zip_files()[0]


    if Path(zip_path).exists():
        print(f"The file {zip_path} exists.")
    else:
        print(f"The file {zip_path} does not exist.")
    # create a temp directory

    tempdir = Path(tempfile.mkdtemp())
    run_string_gsea_bfabric(zip_path, workunit_id, fdr = fdr, base_dir = tempdir)
