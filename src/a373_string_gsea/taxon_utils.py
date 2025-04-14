import polars as pl
import zipfile
from pathlib import Path
from loguru import logger

# Define paths relative to the project root
BASE_DIR = Path(__file__).resolve().parent.parent.parent  # moves two levels up
DATA_MAPPINGS_DIR = BASE_DIR / 'data' / 'mappings'
logger.info(f"Data mappings directory: {DATA_MAPPINGS_DIR}")
SPECIES_ZIP_PATH = DATA_MAPPINGS_DIR / "species.v12.0.zip"
SPECIES_FILE_NAME = "species.v12.0.txt"
NCBI_ZIP_PATH = DATA_MAPPINGS_DIR / "NCBI_nodes.zip"
NCBI_FILE_NAME = "nodes.tsv" # Assuming this is the correct file name inside NCBI_nodes.zip

def read_species_string_data(zip_path: Path = SPECIES_ZIP_PATH, filename: str = SPECIES_FILE_NAME) -> pl.DataFrame | None:
    """Reads the species data from the STRING zip file.
    Args:
        zip_path (Path): The path to the zip file.
        filename (str): The name of the file within the zip archive.

    Returns:
        pl.DataFrame: A Polars DataFrame containing the species data.
    """
    with zipfile.ZipFile(zip_path, 'r') as zip_file:
        # Open and read the file content using Polars.
        with zip_file.open(filename) as file_obj:
            df = pl.read_csv(file_obj, encoding="utf8", truncate_ragged_lines=True, separator="\t")
            df = df.rename({"#taxon_id": "taxon_id"})
            return df


def read_ncbi_nodes_data(zip_path: Path = NCBI_ZIP_PATH, filename: str = NCBI_FILE_NAME) -> pl.DataFrame | None:
    with zipfile.ZipFile(zip_path, 'r') as zip_file:
        # Open and read the file content using Polars.
        with zip_file.open(filename) as file_obj:
            df = pl.read_csv(file_obj, encoding="utf8", separator="\t")
            return df


def get_organism_for_string(spec_id: int, ncbi: pl.DataFrame, species_string: pl.DataFrame) -> int | None:
    """
    Given an initial taxon ID, recursively find an identifier that is present in the species_string
    DataFrame. If the spec_id is not directly present in species_string's 'X.taxon_id' column, the function
    retrieves the corresponding row from ncbi to obtain the parent_taxon_id. It then checks if that parent's ID
    is in species_string. If not, the process continues recursively.

    Args:
        spec_id (int): The initial species taxon ID.
        ncbi (pl.DataFrame): A Polars DataFrame containing at least the columns 'taxon_id' and 'parent_taxon_id'.
        species_string (pl.DataFrame): A Polars DataFrame containing the species IDs in the column 'X.taxon_id'.

    Returns:
        int | None: The species ID found in species_string, or None if no suitable ancestor is found.
    """
    # Check if the current spec_id is in species_string.
    if spec_id in species_string["taxon_id"].to_list():
        return spec_id

    # Filter the ncbi dataframe to get the row for the given spec_id.
    get_taxon = ncbi.filter(pl.col("taxon_id") == spec_id)
    if get_taxon.height == 0:
        logger.error(f"Taxon id {spec_id} not found in ncbi dataframe.")
        return None

    # Extract the parent's taxon id from the returned row.
    parent_taxon_id = get_taxon["parent_taxon_id"].to_list()[0]

    # Check if the parent's taxon id is in species_string.
    if parent_taxon_id in species_string["taxon_id"].to_list():
        return parent_taxon_id
    else:
        logger.info(f"{parent_taxon_id} is not in String")
        if parent_taxon_id == 1:
            return None
        else:
            # Recursively call the function with the parent's taxon id.
            return get_organism_for_string(parent_taxon_id, ncbi, species_string)


# Example Usage (can be removed or kept for testing)
if __name__ == "__main__":
    logger.info("Running example usage...")
    
    # Load data
    string_species = read_species_string_data()
    ncbi_nodes = read_ncbi_nodes_data()

    if string_species is not None and ncbi_nodes is not None:
        logger.info("Data loaded successfully.")
        
        # Test cases
        test_ids = [
            559292, # Yeast - should be found directly
            4932,   # Saccharomyces cerevisiae (parent of 559292) - should be found directly if present, or via parent logic if 559292 is used
            9606,   # Human - should be found directly
            10090,  # Mouse - should be found directly
            511145, # E. coli K12 - should be found directly
            83333,  # E. coli (generic parent) - check if parent (e.g., 511145) is valid
            123456789 # Non-existent ID
        ]

        for test_id in test_ids:
            logger.info(f"\n--- Testing Taxon ID: {test_id} ---")
            valid_id = get_organism_for_string(test_id, ncbi_nodes, string_species)
            if valid_id is not None:
                logger.success(f"Result for {test_id}: Found valid STRING ID -> {valid_id}")
            else:
                logger.warning(f"Result for {test_id}: No valid STRING ID found (checked direct and parent).")
    else:
        logger.error("Failed to load necessary data for example usage.")
