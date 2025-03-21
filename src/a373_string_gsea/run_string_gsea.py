#%% md

import re
import zipfile
import io
import glob
import yaml
from loguru import logger
from a373_string_gsea.stringgsea import StringGSEA
from collections import Counter
from pathlib import Path
import subprocess



if False:
    URL = "https://version-12-0.string-db.org/api/json/get_api_key"
    response = requests.get(URL)
    print("Status code:", response.status_code)
    print("Returned data:", response.json()['api_key'])
# create a test for this function please

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

def get_oxes(zip_path : str):
    oxes = {}
    with zipfile.ZipFile(zip_path, 'r') as z:
        fasta_files = [f for f in z.namelist() if f.endswith(('.fas', '.fasta'))]
        for file in fasta_files:
            with z.open(file) as f:
                file_bytes = f.read()
                ox = get_ox_fields(io.BytesIO(file_bytes))
                oxes[file] = ox
    return oxes

def find_zip_files():
    pattern = re.compile(r'^(\d{7}|DEA_[^/]+)\.zip$')
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

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result


if __name__ == '__main__':
    species: int = 9606
    api_key = "b36F8oaRJwFZ"
    workunit_id = "322025"


    workunit_id =  str(extract_workunit_id_from_file("params.yml"))
    logger.info(f"Workunit ID: {workunit_id}")
    fdr: float = 0.25

    # Path to your zip file
    # zip_path = Path(__file__).parents[1]/'2806039.zip'
    zip_path = find_zip_files()[0]
    res = get_oxes(zip_path)
    merged = [item for sublist in res.values() for item in sublist]
    counter = Counter(merged)
    # Get the most common element (returns a list of tuples)
    species = counter.most_common(1)[0][0]

    dataframes = StringGSEA.get_rank_files(zip_path)
    gsea = StringGSEA(api_key, workunit_id, dataframes, species, fdr)

    gsea.string_gsea()
    logger.info(f"Job submitted successfully.{gsea.res_job_id}")

    #gsea.res_job_id = {'C37638WU322006/Bait_FAN1~FAN1.rnk': 'bhDZz6P6jgmx'}
    gsea.pull_results()
    logger.info("got results")
    path = gsea.zip_folder(gsea.write_results())
    outputs_yml(path)
    logger.info("zipped results.")
    status = gsea.save_link()
    logger.info("saved link")
    



