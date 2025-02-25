import shlex
import zipfile
import io
import polars as pl
import json
import requests
import time
import subprocess
from loguru import logger
from pathlib import Path
import shutil


class StringGSEA:

    def __init__(self, api_key:str,
                 workunit_id:str,
                 dataframes:dict,
                 species:int = 9606,
                 fdr:float=0.25):
        self.api_key =  api_key
        self.workunit_id = workunit_id
        self.species = species
        self.fdr = fdr
        self.dataframes = dataframes

    @staticmethod
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

    def _string_gsea(self, rank_data:str) -> str:
        string_api_url = "https://version-12-0.string-db.org/api"
        output_format = "json"
        method = "valuesranks_enrichment_submit"
        request_url = "/".join([string_api_url, output_format, method])
        params = {
            "species": self.species,  # NCBI/STRING species identifier (e.g., 9606 for human)
            "caller_identity": "www.fgcz.ch",
            "identifiers": rank_data,
            "api_key": self.api_key,
            "ge_fdr": self.fdr,
            "ge_enrichment_rank_direction": -1
        }
        response = requests.post(request_url, data=params)
        data = json.loads(response.text)[0]
        if 'status' in data and data['status'] == 'error':
            logger.error("Status:", data['status'])
            logger.error("Message:", data['message'])
        else:
            job_id = data["job_id"]
            logger.info(f"Job submitted successfully. Job ID: {job_id}")
            return job_id

    def string_gsea(self):
        # submits gsea job to string-db
        self.res_job_id = {}
        for name, out in self.dataframes.items():
            output_str = out.write_csv(separator="\t", include_header=False)
            self.res_job_id[name] = self._string_gsea(output_str)
        return self


    def _pull_results(self, job_id:str, sleep_t:int = 10) -> dict:
        URL = f"https://version-12-0.string-db.org/api/json/valuesranks_enrichment_status?api_key={self.api_key}&job_id={job_id}"
        while True:
            response = requests.get(URL)
            data = response.json()
            status = data[0].get('status')
            logger.info(f"Status: {status}")
            if status == "finished":
                break
            time.sleep(sleep_t)
        return data[0]

    def pull_results(self):
        # queries if resdata is completed.
        self.res_data = {}
        for name, job_id in self.res_job_id.items():
            logger.info(f"job_id: {job_id}")
            logger.info(f"api_key: {self.api_key}")
            self.res_data[name] = self._pull_results(job_id)
        return self

    def _save_link(self, link:str, name:str) -> dict:
        # Define your arguments
        cmd = ["bfabric_save_link_to_workunit.py", self.workunit_id, link, name]
        logger.info(shlex.join(cmd))
        # Run the command and capture the output
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the command executed successfully
        if result.returncode == 0:
            return {result.returncode : result.stdout}
        else:
            return {result.returncode : result.stderr}

    def save_link(self) -> dict:
        save_status = {}
        for name, data in self.res_data.items():
            link_url = data['page_url']
            p = Path(name)
            link_name = "Result GSEA String for Bait [ " + p.stem + "]"
            save_status[name] = self._save_link(link_url, link_name)
        return save_status

    def write_results(self, path:Path = Path(".")):
        res_path = path / f"WU_{self.workunit_id}_GSEA"
        res_path.mkdir(exist_ok=True)

        for name, data in self.res_data.items():
            path = Path(name)
            logger.info(f"Results for {name}:")
            download_url = data['download_url']
            response = requests.get(download_url)
            response.raise_for_status()  # Raises an error for bad responses
            # Option 1: Save the TSV content to a file
            with open(res_path / path.with_suffix(".tsv").name, 'wb') as f:
                f.write(response.content)

            image_url = data['graph_url']
            response = requests.get(image_url)
            response.raise_for_status()  # Raises an error for bad responses
            with open(res_path / path.with_suffix(".png").name, 'wb') as f:
                f.write(response.content)

        return res_path

    @staticmethod
    def zip_folder(folder_path: Path) -> Path:
        folder_path = folder_path.resolve()
        archive_base = folder_path.parent / folder_path.name
        zip_file_path = shutil.make_archive(str(archive_base), 'zip', root_dir=str(folder_path))
        return Path(zip_file_path)

