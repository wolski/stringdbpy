import json
import os
import requests
import time
from loguru import logger
from pathlib import Path
import shutil

from a373_string_gsea.gsea_utilities import get_rank_files
from a373_string_gsea.get_species import OxFieldsZip

class StringGSEA:

    def __init__(self, api_key:str,
                 workunit_id:str,
                 rank_dataframes:dict,
                 species:int = 9606,
                 fdr:float=0.25):
        self.api_key =  api_key
        self.workunit_id = workunit_id
        self.species = species
        self.fdr = fdr
        self.rank_dataframes = rank_dataframes
        self.res_data = {}
        self.res_job_id = {}

    @classmethod
    def from_serialized_json(cls, json_path: Path):
        """Initializes a StringGSEA instance from a serialized JSON session file.

        Reads configuration and results data from the JSON.

        Args:
            json_path: Path to the serialized JSON file ('serialized_gsea_session.json').

        Returns:
            A StringGSEA instance populated with data from the JSON file.

        Raises:
            FileNotFoundError: If the json_path does not exist.
            json.JSONDecodeError: If the file is not valid JSON.
            KeyError: If the JSON structure is missing 'config' or 'results', or keys within them.
            ValueError: If a string key in results cannot be split into a tuple key.
        """
        if not json_path.is_file():
            raise FileNotFoundError(f"JSON file not found at {json_path}")

        try:
            with open(json_path, 'r') as f:
                session_data = json.load(f)
        except json.JSONDecodeError as e:
            logger.error(f"Failed to decode JSON from {json_path}: {e}")
            raise

        # Extract config
        try:
            config = session_data['config']
            api_key = config['api_key']
            workunit_id = config['workunit_id']
            species = config['species']
            fdr = config['fdr']
        except KeyError as e:
            logger.error(f"Missing expected configuration key in {json_path}: {e}")
            raise

        # Create a new instance using the extracted config
        # Initialize dataframes as empty dict, as it's not stored
        instance = cls(api_key=api_key,
                       workunit_id=workunit_id,
                       rank_dataframes={}, # dataframes are not stored/restored
                       species=species,
                       fdr=fdr)

        # Reconstruct res_data with tuple keys
        reconstructed_res_data = {}
        try:
            serialized_results = session_data['results']
            for str_key, value in serialized_results.items():
                try:
                    # Split by the "~" to handle potential underscores in keys
                    parts = str_key.rsplit('~', 1)
                    if len(parts) == 2:
                        outer_key, inner_key = parts
                        tuple_key = (outer_key, inner_key)
                        reconstructed_res_data[tuple_key] = value
                    else:
                        # Handle cases where splitting might fail
                        raise ValueError(f"Could not reconstruct tuple key from string key: '{str_key}'")
                except Exception as e:
                    logger.error(f"Error processing results key '{str_key}': {e}")
                    raise # Re-raise after logging
        except KeyError:
             logger.error(f"Missing 'results' key in {json_path}")
             raise


        instance.res_data = reconstructed_res_data
        # Note: res_job_id is not stored/restored from the JSON.
        # This instance is primarily for accessing results.

        logger.info(f"StringGSEA instance initialized from {json_path}")
        return instance

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

        if response.ok :
            data = json.loads(response.text)[0]
            if 'status' in data and data['status'] == 'error':
                logger.error("Status:", data['status'])
                logger.error("Message:", data['message'])
            else:
                job_id = data["job_id"]
                logger.info(f"Job submitted successfully. Job ID: {job_id}")
                return job_id
        else:
            logger.info(response.text)
            logger.error(f"Request failed with status code: {response.status_code}")
            raise RuntimeError("Request failed with status code: {}".format(response.status_code))

    def submit(self):
        # submits gsea job to string-db
        for name, out in self.rank_dataframes.items():
            output_str = out.write_csv(separator="\t", include_header=False)
            self.res_job_id[name] = self._string_gsea(output_str)
        return self


    def _pull_results(self, job_id:str, sleep_t:int = 10, max_time = 3600) -> dict:
        self.URL = f"https://version-12-0.string-db.org/api/json/valuesranks_enrichment_status?api_key={self.api_key}&job_id={job_id}"
        elapsed = 0
        data = None
        while elapsed < max_time:
            response = requests.get(self.URL)
            data = response.json()
            status = data[0].get('status')
            logger.info(f"time : {elapsed} Status: {status}")
            if status == "success":
                break
            if status == "nothing found":
                 # Raise an exception if no results were found.
                raise Exception(f"No results found for job_id {job_id}")
            if status == "unknown organism":
                # Raise an exception if the organism is unknown.
                raise Exception(f"Unknown organism encountered for job_id {job_id}")
            time.sleep(sleep_t)
            elapsed += sleep_t
        return data[0]

    def pull_results(self):
        for name, job_id in self.res_job_id.items():
            logger.info(f"job_id: {job_id}")
            logger.info(f"api_key: {self.api_key}")
            self.res_data[name] = self._pull_results(job_id)
        return self

    def get_res_path(self, path : Path = Path(".")) -> Path:
        res_path = path / f"WU_{self.workunit_id}_GSEA"
        res_path.mkdir(exist_ok=True)
        return res_path

    def get_links(self) -> dict:
        links_dict = {}
        for (outer_key, inner_key), data in self.res_data.items():
            if data.get('status') == "success" and 'page_url' in data:
                # Initialize the inner dictionary if it doesn't exist
                if outer_key not in links_dict:
                    links_dict[outer_key] = {}
                
                # Add the link to the inner dictionary
                links_dict[outer_key][inner_key] = data['page_url']
                logger.info(f"Added link for {outer_key}/{inner_key}")
        
        return links_dict
    
    def _write_links_to_file(self, outer_key: str, path: Path = Path(".")) -> Path:
        links_dict = self.get_links()
        
        if outer_key not in links_dict:
            logger.warning(f"No links found for outer_key: {outer_key}")
            return None
        
        subfolder_path = path / str(outer_key)
        subfolder_path.mkdir(exist_ok=True)
        
        filename = "links.txt"
        file_path = subfolder_path / filename
        
        with open(file_path, 'w') as f:
            for inner_key, link_url in links_dict[outer_key].items():
                f.write(f"{inner_key}: {link_url}\n")
        
        logger.info(f"Wrote links for {outer_key} to {file_path}")
        return file_path


    def write_links(self, path: Path = Path(".")) -> dict:
        path = self.get_res_path(path)
        written_files = {}
        links_dict = self.get_links()
        
        # Iterate over each outer_key in the links dictionary.
        for outer_key in links_dict:
            file_path = self._write_links_to_file(outer_key, path)
            if file_path is not None:
                written_files[outer_key] = file_path
            else:
                logger.warning(f"No file was written for outer_key: {outer_key}")
        return written_files


    def write_rank_files(self, path: Path = Path(".")):
        res_path = self.get_res_path(path)
        res_path.mkdir(exist_ok=True)

        for (outer_key, inner_key), df in self.rank_dataframes.items():
            subfolder_path = res_path / str(outer_key)
            subfolder_path.mkdir(exist_ok=True)

            filename = f"{inner_key}.rnk"
            logger.info(f"Writing rank file for {outer_key}/{inner_key} -> {filename}")

            file_path = subfolder_path / filename
            with open(file_path, 'wb') as f:
                f.write(df.write_csv(separator="\t", include_header=False).encode())
        return res_path


    def write_gsea_results(self, path:Path = Path(".")) :
        res_path = self.get_res_path(path)

        for (outer_key, inner_key), data in self.res_data.items():
            if data.get('status') == "success":
                logger.info(f"Results for {outer_key} {inner_key}:")
                download_url = data['download_url']
                response = requests.get(download_url)
                response.raise_for_status()  # Raises an error for bad responses
                # Option 1: Save the TSV content to a file

                # write downloaded tsv dast
                subfolder_path = res_path / str(outer_key)
                subfolder_path.mkdir(exist_ok=True)
                filename = f"{inner_key}.tsv"
                file_path = subfolder_path / filename
                with open(file_path, 'wb') as f:
                    f.write(response.content)

                # Write downloaded image
                filename = f"{inner_key}.png"
                file_path = subfolder_path / filename
                image_url = data['graph_url']
                response = requests.get(image_url)
                response.raise_for_status()  # Raises an error for bad responses
                with open(file_path, 'wb') as f:
                    f.write(response.content)
        return res_path

    def serialize_results(self, path: Path = Path(".")) -> Path:
        """Serializes the self.res_data dictionary and configuration to a JSON file.

        Converts tuple keys in res_data to strings before serialization.
        Stores configuration (api_key, workunit_id, species, fdr) as well.

        Args:
            path: The base directory to save the file in. Defaults to current directory.

        Returns:
            The Path object of the created JSON file.
        """
        res_path = self.get_res_path(path)
        filename = "serialized_gsea_session.json" # Renamed for clarity
        file_path = res_path / filename

        # Convert tuple keys in res_data to strings for JSON compatibility
        serializable_results = {}
        for key_tuple, value in self.res_data.items():
            # Assuming key_tuple is like (outer_key, inner_key)
            str_key = f"{key_tuple[0]}~{key_tuple[1]}"
            serializable_results[str_key] = value

        # Combine config and results
        output_data = {
            "config": {
                "api_key": self.api_key, # Consider security implications if sensitive
                "workunit_id": self.workunit_id,
                "species": self.species,
                "fdr": self.fdr
            },
            "results": serializable_results
        }

        try:
            with open(file_path, 'w') as f:
                json.dump(output_data, f, indent=4)
            logger.info(f"Serialized GSEA session written to {file_path}")
            return file_path
        except TypeError as e:
            logger.error(f"Failed to serialize results data to JSON: {e}")
            raise
        except IOError as e:
            logger.error(f"Failed to write serialized results to file {file_path}: {e}")
            raise
        return file_path


    @staticmethod
    def zip_folder(folder_path: Path) -> Path:
        folder_path = folder_path.resolve()
        archive_base = folder_path.parent / folder_path.name
        zip_file_path = shutil.make_archive(str(archive_base), 'zip', root_dir=str(folder_path))
        return Path(zip_file_path)


if __name__ == '__main__':
    api_key = "b36F8oaRJwFZ"
    fdr: float = 0.25
    current_directory = os.getcwd()
    print("Current working directory:", current_directory)

    # zip_path = "./testing_examples/fromDEA/2790797.zip"
    zip_path = "./tests/data/DE_mouse_fasta_rnk.zip"
    workunit_id = "abcd"
    species = OxFieldsZip.get_species_from_oxes(zip_path)
    logger.info(f"species is : {species}")
    dataframes = get_rank_files(zip_path)

    gsea = StringGSEA(api_key, workunit_id, dataframes, species, fdr)
    gsea.submit()
    logger.info(f"Job submitted successfully.{gsea.res_job_id}")
    gsea.pull_results()
    logger.info("got results")
    p = Path("./tests/data")
    p.exists()
    serres_file = gsea.serialize_results(p)
    
    gsea.write_rank_files()
    gsea.write_links()
    result_dir = gsea.write_gsea_results()
    serres_file = "/Users/witoldwolski/__checkout/slurmworker/config/A373_STRING_GSEA/WU_abcd_GSEA/serialized_gsea_session.json"
    gsea2 = StringGSEA.from_serialized_json(serres_file)
    gsea2.workunit_id = "xyz"
    gsea.write_gsea_results()