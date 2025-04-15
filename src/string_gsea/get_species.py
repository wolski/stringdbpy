import io
import re # Added missing import
from collections import Counter # Corrected import
import requests
from loguru import logger
import polars as pl
import zipfile
from pathlib import Path
from string_gsea.gsea_utilities import get_rank_files


class GetTaxonID:
    @staticmethod
    def _fetch_ncbi_taxon_id(identifier):
        """
        Fetches the ncbiTaxonId for a given identifier from the STRING API.
        
        Parameters:
            identifier (str): The identifier to query (e.g., "P05737").
        
        Returns:
            int or None: The ncbiTaxonId if found, otherwise None.
        """
        url = "https://string-db.org/api/json/get_string_ids"
        params = {"identifiers": identifier}
        
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        
        # Check if the response contains a list and return the "ncbiTaxonId" from the first element.
        if data and isinstance(data, list):
            return data[0].get("ncbiTaxonId")
        return None

    @staticmethod
    def determine_species(df: pl.DataFrame, nr: int = 10) -> int:
        """
        Determines the most likely species from a rank file DataFrame by sampling identifiers
        and querying the STRING API for their NCBI Taxon IDs.

        Parameters:
            df (pl.DataFrame): DataFrame read from a rank file. Assumes the first column contains identifiers.
            nr (int): The number of identifiers to sample. Defaults to 100.

        Returns:
            str: The most frequent ncbiTaxonId found.

        Raises:
            ValueError: If the DataFrame is empty, has no valid identifiers,
                        or if no taxon IDs could be retrieved.
        """

        # Check that the input DataFrame is valid.
        if df.is_empty() or df.shape[1] == 0:
            raise ValueError("Input DataFrame is empty or has no columns.")

        # Get identifiers from the first column, removing nulls and duplicates.
        identifier_col_name = df.columns[0]
        identifiers = df.select(pl.col(identifier_col_name)).drop_nulls().unique()
        if identifiers.height == 0:
            raise ValueError("No valid identifiers found in the first column.")

        # Determine the sample size and sample identifiers.
        sample_size = min(nr, identifiers.height)
        sampled_ids = identifiers.sample(n=sample_size, shuffle=True).to_series().to_list()

        # Fetch taxon IDs for sampled identifiers.
        taxon_ids = []
        for identifier in sampled_ids:
            try:
                taxon_id = GetTaxonID._fetch_ncbi_taxon_id(identifier)
                if taxon_id is not None:
                    taxon_ids.append(taxon_id)
            except Exception:
                # Skip identifiers that raise an exception during API call.
                continue

        # If no taxon IDs were fetched, raise an exception.
        if not taxon_ids:
            raise ValueError("Could not fetch any taxon IDs for the sampled identifiers.")

        # Count frequencies and return the most common taxon ID.
        most_common_taxon = Counter(taxon_ids).most_common(1)[0][0]
        return most_common_taxon


class OxFieldsZip:
    @staticmethod
    def get_ox_fields(fasta_content: io.BytesIO) -> list[str]:
        pattern = re.compile(r'OX=(\d+)') # re was not defined before
        ox_values = []

        fasta_content = io.TextIOWrapper(fasta_content, encoding='utf-8')
        for line in fasta_content:
            # Check only header lines starting with '>'
            if line.startswith('>'):
                match = pattern.search(line)
                if match:
                    ox_values.append(match.group(1))

        return ox_values
    
    @staticmethod
    def get_species_from_oxes(zip_path : str) -> int:
        oxes = {}
        with zipfile.ZipFile(zip_path, 'r') as z:
            fasta_files = [f for f in z.namelist() if f.endswith(('.fas', '.fasta'))]
            for file in fasta_files:
                with z.open(file) as f:
                    file_bytes = f.read()
                    ox = OxFieldsZip.get_ox_fields(io.BytesIO(file_bytes))
                    oxes[file] = ox

        merged = [item for sublist in oxes.values() for item in sublist]
        counter = Counter(merged)
        # Get the most common element (returns a list of tuples)
        species = counter.most_common(1)[0][0]

        return int(species)




class TaxonUtils:
    """A class for handling taxon-related utilities and data processing."""
    
    def __init__(self):
        # Define paths relative to the project root
        self.base_dir = Path(__file__).resolve().parent.parent.parent  # moves two levels up
        self.data_mappings_dir = self.base_dir / 'data' / 'mappings'
        logger.info(f"Data mappings directory: {self.data_mappings_dir}")
        self.species_zip_path = self.data_mappings_dir / "species.v12.0.zip"
        self.species_file_name = "species.v12.0.txt"
        self.ncbi_zip_path = self.data_mappings_dir / "NCBI_nodes.zip"
        self.ncbi_file_name = "nodes.tsv"  # Assuming this is the correct file name inside NCBI_nodes.zip
        
        # Load data during initialization
        self.string_species_df = self._read_species_string_data()
        self.ncbi_nodes_df = self._read_ncbi_nodes_data()
        
        if self.string_species_df is None or self.ncbi_nodes_df is None:
            logger.error("Failed to load necessary data during initialization.")
    
    def _read_species_string_data(self) -> pl.DataFrame | None:
        """Reads the species data from the STRING zip file.

        Returns:
            pl.DataFrame: A Polars DataFrame containing the species data.
        """
        with zipfile.ZipFile(self.species_zip_path, 'r') as zip_file:
            # Open and read the file content using Polars.
            with zip_file.open(self.species_file_name) as file_obj:
                df = pl.read_csv(file_obj, encoding="utf8", truncate_ragged_lines=True, separator="\t")
                df = df.rename({"#taxon_id": "taxon_id"})
                return df

    def _read_ncbi_nodes_data(self) -> pl.DataFrame | None:
        """Reads the NCBI nodes data from the zip file.

        Returns:
            pl.DataFrame: A Polars DataFrame containing the NCBI nodes data.
        """
        with zipfile.ZipFile(self.ncbi_zip_path, 'r') as zip_file:
            # Open and read the file content using Polars.
            with zip_file.open(self.ncbi_file_name) as file_obj:
                df = pl.read_csv(file_obj, encoding="utf8", separator="\t")
                return df

    def get_organism_for_string(self, spec_id: int) -> int | None:
        """
        Given an initial taxon ID, recursively find an identifier that is present in the species_string
        DataFrame. If the spec_id is not directly present in species_string's 'X.taxon_id' column, the function
        retrieves the corresponding row from ncbi to obtain the parent_taxon_id. It then checks if that parent's ID
        is in species_string. If not, the process continues recursively.

        Args:
            spec_id (int): The initial species taxon ID.

        Returns:
            int | None: The species ID found in species_string, or None if no suitable ancestor is found.
        """
        # Check if the current spec_id is in species_string.
        if spec_id in self.string_species_df["taxon_id"].to_list():
            return spec_id

        # Filter the ncbi dataframe to get the row for the given spec_id.
        get_taxon = self.ncbi_nodes_df.filter(pl.col("taxon_id") == spec_id)
        if get_taxon.height == 0:
            logger.error(f"Taxon id {spec_id} not found in ncbi dataframe.")
            return None

        # Extract the parent's taxon id from the returned row.
        parent_taxon_id = get_taxon["parent_taxon_id"].to_list()[0]

        # Check if the parent's taxon id is in species_string.
        if parent_taxon_id in self.string_species_df["taxon_id"].to_list():
            return parent_taxon_id
        else:
            logger.info(f"{parent_taxon_id} is not in String")
            if parent_taxon_id == 1:
                return None
            else:
                # Recursively call the function with the parent's taxon id.
                return self.get_organism_for_string(parent_taxon_id)


def get_species_taxon(zip_path : Path, nr = 10) -> int:
    taxon = OxFieldsZip.get_species_from_oxes(str(zip_path))
    taxon2 = TaxonUtils().get_organism_for_string(taxon)
    if taxon2 is not None:
        logger.info(f"Taxon : {taxon2} found in species_string")
        return taxon2
    else:
        df_list = get_rank_files(YEAST_RNK_ZIP_PATH)
        df = list(df_list.values())[0]

        taxon3 = GetTaxonID.determine_species(df, nr = nr)
        logger.info(f"{taxon3} retrieved from string by searching {nr} identifiers.")

        return taxon3


# Example Usage (can be removed or kept for testing)
if __name__ == "__main__":
    TEST_DATA_DIR = Path(__file__).resolve().parent.parent.parent / 'tests/data'
    YEAST_RNK_ZIP_PATH = TEST_DATA_DIR / 'DE_yeast_fasta_rnk.zip' # Path for yeast rank file zip

    p = Path(YEAST_RNK_ZIP_PATH)
    taxon = get_species_taxon(p, nr=5)
    logger.info(f"Taxon ID: {taxon}")

    