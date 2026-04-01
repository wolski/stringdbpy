"""Species mapping via embedded STRING-DB and NCBI taxonomy data."""

import importlib.resources
import zipfile

import polars as pl
from loguru import logger


class TaxonUtils:
    """Resolves taxon IDs to STRING-DB supported species via NCBI taxonomy tree traversal."""

    def __init__(self):
        self.mappings_package = "string_gsea.data.mappings"
        self.species_zip_name = "species.v12.0.zip"
        self.species_file_name = "species.v12.0.txt"
        self.ncbi_zip_name = "NCBI_nodes.zip"
        self.ncbi_file_name = "nodes.tsv"

        self.string_species_df = self._read_species_string_data()
        self.ncbi_nodes_df = self._read_ncbi_nodes_data()

    def _read_species_string_data(self) -> pl.DataFrame:
        """Reads the species data from the STRING zip file."""
        with importlib.resources.path(self.mappings_package, self.species_zip_name) as zip_path:
            with zipfile.ZipFile(zip_path, "r") as zip_file:
                with zip_file.open(self.species_file_name) as file_obj:
                    df = pl.read_csv(
                        file_obj,
                        encoding="utf8",
                        truncate_ragged_lines=True,
                        separator="\t",
                    )
                    df = df.rename({"#taxon_id": "taxon_id"})
                    return df

    def _read_ncbi_nodes_data(self) -> pl.DataFrame:
        """Reads the NCBI nodes data from the zip file."""
        with importlib.resources.path(self.mappings_package, self.ncbi_zip_name) as zip_path:
            with zipfile.ZipFile(zip_path, "r") as zip_file:
                with zip_file.open(self.ncbi_file_name) as file_obj:
                    df = pl.read_csv(file_obj, encoding="utf8", separator="\t")
                    return df

    def get_organism_for_string(self, spec_id: int) -> int | None:
        """
        Given a taxon ID, walk the NCBI taxonomy tree upward until a STRING-DB
        supported species is found.

        Args:
            spec_id: The initial species taxon ID.

        Returns:
            The STRING-DB species ID, or None if no supported ancestor is found.
        """
        if spec_id in self.string_species_df["taxon_id"].to_list():
            return spec_id

        get_taxon = self.ncbi_nodes_df.filter(pl.col("taxon_id") == spec_id)
        if get_taxon.height == 0:
            logger.error(f"Taxon id {spec_id} not found in ncbi dataframe.")
            return None

        parent_taxon_id = get_taxon["parent_taxon_id"].to_list()[0]

        if parent_taxon_id in self.string_species_df["taxon_id"].to_list():
            return parent_taxon_id
        else:
            logger.info(f"{parent_taxon_id} is not in String")
            if parent_taxon_id == 1:
                return None
            else:
                return self.get_organism_for_string(parent_taxon_id)
