#%% md

import re
import zipfile
import io
import glob
import shlex

import yaml
from loguru import logger

from a373_string_gsea import postprocess
from a373_string_gsea.stringgsea import StringGSEA
from collections import Counter
from pathlib import Path
import subprocess
import polars as pl

def get_rank_files(zip_path):
    dataframes = {}
    with zipfile.ZipFile(zip_path, 'r') as z:
        # Filter out the .rnk files from the archive
        rnk_files = [f for f in z.namelist() if f.endswith('.rnk')]
        for file in rnk_files:
            with z.open(file) as f:
                # Read the file content into memory
                file_bytes = f.read()
                # Wrap bytes in a BytesIO stream for polars to read from
                # Adjust the separator (sep) if your file is not tab-delimited
                df = pl.read_csv(io.BytesIO(file_bytes), separator="\t", has_header=False)
                dataframes[file] = df
    return dataframes


def get_ox_fields(fasta_content):
    pattern = re.compile(r'OX=(\d+)')
    ox_values = []

    if isinstance(fasta_content, (bytes, bytearray)):
        fasta_content = io.BytesIO(fasta_content)

    if isinstance(fasta_content, io.BytesIO):
        fasta_content = io.TextIOWrapper(fasta_content, encoding='utf-8')

    for line in fasta_content:
        # Check only header lines starting with '>'
        if line.startswith('>'):
            match = pattern.search(line)
            if match:
                ox_values.append(match.group(1))

    return ox_values

def get_species_from_oxes(zip_path : str):
    oxes = {}
    with zipfile.ZipFile(zip_path, 'r') as z:
        fasta_files = [f for f in z.namelist() if f.endswith(('.fas', '.fasta'))]
        for file in fasta_files:
            with z.open(file) as f:
                file_bytes = f.read()
                ox = get_ox_fields(io.BytesIO(file_bytes))
                oxes[file] = ox

    merged = [item for sublist in oxes.values() for item in sublist]
    counter = Counter(merged)
    # Get the most common element (returns a list of tuples)
    species = counter.most_common(1)[0][0]

    return species

def find_zip_files():
    pattern = re.compile(r'^(\d{7}.*|DEA_[^/]+)\.zip$')
    files = [f for f in glob.glob("*.zip") if pattern.match(f)]
    return files

def extract_workunit_id_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            data = yaml.safe_load(file)
        # Navigate to the 'registration' block and extract 'workunit_id'
        return data.get('registration', {}).get('workunit_id')
    except Exception as e:
        logger.error(f"Error reading or parsing the YAML file: {e}")
        return None

def outputs_yml(search_zip : Path, outputs_yml = "outputs.yml"):
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
        
        with open(outputs_yml, 'w') as file:
            yaml.dump(data, file, default_flow_style=False)
        logger.info(f"YAML file {outputs_yml} has been generated.")


def _save_link(link:str, name:str, workunit_id: str) -> dict:
    # Define your arguments
    cmd = ["bfabric_save_link_to_workunit.py", workunit_id, link, name]
    logger.info(shlex.join(cmd))
    # Run the command and capture the output
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check if the command executed successfully
    if result.returncode == 0:
        return {result.returncode : result.stdout}
    else:
        return {result.returncode : result.stderr}

def save_link(res_data : dict, workunit_id: str) -> dict:
    save_status = {}
    for name, data in res_data.items():
        if data.get('status') == "success":
            link_url = data['page_url']
            p = Path(name)
            link_name = "Result string-db GSEA for [ " + p.stem + "]"
            save_status[name] = _save_link(link_url, link_name,workunit_id)
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


def run_string_gsea(zip_path : str, workunit_id: str, api_key: str, fdr: float = 0.25):

    species = get_species_from_oxes(zip_path)
    dataframes = get_rank_files(zip_path)
    gsea = StringGSEA(api_key, workunit_id, dataframes, species, fdr)
    gsea.string_gsea()
    logger.info(f"Job submitted successfully.{gsea.res_job_id}")
    gsea.pull_results()
    logger.info("got results")
    result_dir = gsea.write_results()
    postprocess.result_to_xlsx(result_dir, workunit_id)
    logger.info(f"Results written to {result_dir}")
    path = gsea.zip_folder(result_dir)
    logger.info(f"Zipped results to {path}")
    outputs_yml(path)
    logger.info("zipped results.")
    status = save_link(gsea.res_data, workunit_id)
    logger.info(f"saved link {status}")


if __name__ == '__main__':
    species: int
    api_key = "b36F8oaRJwFZ"
    workunit_id: str


    workunit_id =  str(extract_workunit_id_from_file("params.yml"))
    logger.info(f"Workunit ID: {workunit_id}")
    fdr: float = 0.25
    zip_path = find_zip_files()[0]
    run_string_gsea(zip_path, workunit_id, api_key, fdr)






