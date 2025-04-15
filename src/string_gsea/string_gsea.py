import json
import os
import requests
import time
from loguru import logger
from pathlib import Path
import shutil
import tempfile

from string_gsea.gsea_utilities import get_rank_files
from string_gsea.get_species import OxFieldsZip

class StringGSEA:

    def __init__(self, api_key:str,
                 workunit_id:str,
                 rank_dataframes:dict,
                 species:int = 9606,
                 fdr:float=0.25,
                 base_path:Path = Path(".")):
        self.api_key =  api_key
        self.workunit_id = workunit_id
        self.species = species
        self.fdr = fdr
        self.rank_dataframes = rank_dataframes
        self.res_data = {}
        self.res_job_id = {}
        self.base_path = base_path

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
            base_path = Path(config.get('base_path', '.'))
        except KeyError as e:
            logger.error(f"Missing expected configuration key in {json_path}: {e}")
            raise

        # Create a new instance using the extracted config
        # Initialize dataframes as empty dict, as it's not stored
        instance = cls(api_key=api_key,
                       workunit_id=workunit_id,
                       rank_dataframes={}, # dataframes are not stored/restored
                       species=species,
                       fdr=fdr,
                       base_path=base_path)

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

    def get_res_path(self) -> Path:
        """Get the path for storing results.
        
        Returns:
            Path: The path for storing results.
        """
        res_path = self.base_path / f"WU_{self.workunit_id}_GSEA"
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
    
    def _write_links_to_file(self, outer_key: str) -> Path:
        """Write links for a specific outer key to a file.
        
        Args:
            outer_key: The outer key for which to write links.
            
        Returns:
            Path: The path to the file where links were written, or None if no links were found.
        """
        links_dict = self.get_links()
        
        if outer_key not in links_dict:
            logger.warning(f"No links found for outer_key: {outer_key}")
            return None
        
        subfolder_path = self.get_res_path() / str(outer_key)
        subfolder_path.mkdir(exist_ok=True)
        
        filename = "links.txt"
        file_path = subfolder_path / filename
        
        with open(file_path, 'w') as f:
            for inner_key, link_url in links_dict[outer_key].items():
                f.write(f"{inner_key}: {link_url}\n")
        
        logger.info(f"Wrote links for {outer_key} to {file_path}")
        return file_path


    def write_links(self) -> dict:
        """Write links for all outer keys to files.
        
        Returns:
            dict: A dictionary mapping outer keys to the paths of the files where links were written.
        """
        written_files = {}
        links_dict = self.get_links()
        
        # Iterate over each outer_key in the links dictionary.
        for outer_key in links_dict:
            file_path = self._write_links_to_file(outer_key)
            if file_path is not None:
                written_files[outer_key] = file_path
            else:
                logger.warning(f"No file was written for outer_key: {outer_key}")
        return written_files


    def write_rank_files(self):
        """Write rank files for all rank dataframes.
        
        Returns:
            Path: The path where rank files were written.
        """
        res_path = self.get_res_path()
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


    def write_gsea_results(self):
        """Write GSEA results for all successful jobs.
        
        Returns:
            Path: The path where GSEA results were written.
        """
        res_path = self.get_res_path()

        for (outer_key, inner_key), data in self.res_data.items():
            if data.get('status') == "success":
                logger.info(f"Results for {outer_key} {inner_key}:")
                download_url = data['download_url']
                response = requests.get(download_url)
                response.raise_for_status()  # Raises an error for bad responses
                
                # Create subfolder for this result
                subfolder_path = res_path / str(outer_key)
                subfolder_path.mkdir(exist_ok=True)
                
                # Save the results
                filename = f"{inner_key}_gsea_results.tsv"
                file_path = subfolder_path / filename
                
                with open(file_path, 'wb') as f:
                    f.write(response.content)
                
                logger.info(f"Saved GSEA results to {file_path}")
        
        return res_path

    def serialize_results(self) -> Path:
        """Serialize the results to a JSON file.
        
        Returns:
            Path: The path to the serialized JSON file.
        """
        # Create a serializable version of the results
        serializable_results = {}
        for (outer_key, inner_key), value in self.res_data.items():
            # Convert tuple keys to strings for JSON serialization
            str_key = f"{outer_key}~{inner_key}"
            serializable_results[str_key] = value
        
        # Create the serialized data structure
        serialized_data = {
            'config': {
                'api_key': self.api_key,
                'workunit_id': self.workunit_id,
                'species': self.species,
                'fdr': self.fdr,
                'base_path': str(self.base_path)
            },
            'results': serializable_results
        }
        
        # Write to file
        output_path = self.base_path / 'serialized_gsea_session.json'
        with open(output_path, 'w') as f:
            json.dump(serialized_data, f, indent=2)
        
        logger.info(f"Serialized GSEA session to {output_path}")
        return output_path


    @staticmethod
    def zip_folder(folder_path: Path) -> Path:
        folder_path = folder_path.resolve()
        archive_base = folder_path.parent / folder_path.name
        zip_file_path = shutil.make_archive(str(archive_base), 'zip', root_dir=str(folder_path))
        return Path(zip_file_path)


if __name__ == '__main__':
    # Define the path to the test data directory relative to the test file
    api_key = "b36F8oaRJwFZ"
    fdr: float = 0.25


    current_directory = os.getcwd()
    print("Current working directory:", current_directory)

    zip_path = "../../tests/data/DE_mouse_fasta_rnk.zip"
    workunit_id = "abcd"
    species = OxFieldsZip.get_species_from_oxes(zip_path)
    logger.info(f"species is : {species}")
    dataframes = get_rank_files(zip_path)

    # Create a temporary directory for testing
    tempdir = Path(tempfile.mkdtemp())
    logger.info(f"Created temporary directory: {tempdir}")

    # Initialize StringGSEA with the temporary directory as base_path
    gsea = StringGSEA(api_key, workunit_id, dataframes, species, fdr, base_path=tempdir)
    gsea.submit()
    logger.info(f"Job submitted successfully.{gsea.res_job_id}")
    gsea.pull_results()
    logger.info("got results")

    # Serialize results
    serialized_json_file = gsea.serialize_results()
    
    # Write files using the base_path set in the constructor
    gsea.write_rank_files()
    gsea.write_links()
    result_dir = gsea.write_gsea_results()

    # Create a new instance from the serialized JSON
    gsea2 = StringGSEA.from_serialized_json(serialized_json_file)

    # Change the workunit_id and write results again
    gsea2.workunit_id = "xyz"
    gsea2.write_gsea_results()
