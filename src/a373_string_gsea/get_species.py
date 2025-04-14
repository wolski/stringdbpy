import io
import re # Added missing import
from collections import Counter # Corrected import
import zipfile
import requests
import polars as pl
from loguru import logger



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
    def determine_species(df: pl.DataFrame, nr: int = 100) -> str:
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
        return str(most_common_taxon)


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
    def get_species_from_oxes(zip_path : str):
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

        return species
