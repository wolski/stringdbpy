"""Low-level HTTP client for STRING-DB GSEA API."""

import time
from collections import Counter
from dataclasses import asdict

import polars as pl
import requests
from loguru import logger

from string_gsea.gsea_config import STRING_API_BASE_DEFAULT, GSEAConfig


class StringAPIClient:
    """Handles HTTP communication with the STRING-DB GSEA API."""

    def __init__(self, config: GSEAConfig, species: int):
        self.config = config
        self.species = species
        self._config_dict = asdict(config)

    def submit_ranks(self, rank_data: str) -> str:
        """Submit rank data to STRING-DB and return the job ID."""
        url = f"{self.config.api_base_url}/json/valuesranks_enrichment_submit"
        params = {
            "species": self.species,
            "caller_identity": self._config_dict["caller_identity"],
            "identifiers": rank_data,
            "api_key": self._config_dict["api_key"],
            "ge_fdr": self._config_dict["fdr"],
            "ge_enrichment_rank_direction": self._config_dict["ge_enrichment_rank_direction"],
        }
        resp = requests.post(url, data=params)
        resp.raise_for_status()
        data = resp.json()[0]
        if data.get("status") == "error":
            msg = data.get("message", "no message")
            logger.error(f"STRING-db error: {msg}")
            raise RuntimeError(f"STRING-db error: {msg}")
        job_id = data["job_id"]
        logger.info(f"Job submitted: {job_id}")
        return job_id

    def poll_job(self, job_id: str, sleep_t: int = 10, max_time: int = 3600) -> dict:
        """Poll a submitted job until success or failure."""
        status_url = f"{self.config.api_base_url}/json/valuesranks_enrichment_status"
        elapsed = 0
        while elapsed < max_time:
            resp = requests.get(status_url, params={"api_key": self._config_dict["api_key"], "job_id": job_id})
            resp.raise_for_status()
            record = resp.json()[0]
            st = record.get("status")
            logger.info(f"Polling {job_id}: status={st} after {elapsed}s")
            if st == "success":
                return record
            if st in ("nothing found", "unknown organism"):
                raise RuntimeError(f"Polling error {st} for job {job_id}")
            time.sleep(sleep_t)
            elapsed += sleep_t
        raise TimeoutError(f"Job {job_id} did not finish in {max_time}s")


# ---------------------------------------------------------------------------
# Module-level functions for species detection via STRING API
# ---------------------------------------------------------------------------


def _fetch_ncbi_taxon_ids(identifiers, api_base_url: str = STRING_API_BASE_DEFAULT):
    """
    Fetch ncbiTaxonId for one or more identifiers from the STRING API.

    Parameters:
        identifiers: A single identifier string or a list of them.
        api_base_url: STRING API base URL.

    Returns:
        dict mapping each input identifier to its ncbiTaxonId (int or None).
    """
    if isinstance(identifiers, (list, tuple)):
        id_string = "\r".join(identifiers)
    else:
        id_string = identifiers

    url = f"{api_base_url}/json/get_string_ids"
    params = {"identifiers": id_string}
    logger.info(f"Fetching ncbiTaxonId for identifiers: {url}?{id_string}")

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    result = {}
    if isinstance(data, list):
        for entry in data:
            input_id = entry.get("queryItem") or entry.get("inputId") or None
            if input_id is None:
                parts = id_string.split("\r")
                input_id = parts[len(result)]
            result[input_id] = entry.get("ncbiTaxonId")
    return result


def _sample_identifiers(df: pl.DataFrame, nr: int = 10) -> list[str]:
    """Sample unique identifiers from the first column of a DataFrame.

    Raises:
        ValueError: If the DataFrame is empty or has no valid identifiers.
    """
    if df.is_empty() or df.shape[1] == 0:
        raise ValueError("Input DataFrame is empty or has no columns.")

    identifier_col = df.columns[0]
    identifiers = df.select(pl.col(identifier_col)).drop_nulls().unique()
    if identifiers.height == 0:
        raise ValueError("No valid identifiers found in the first column.")

    sample_size = min(nr, identifiers.height)
    return identifiers.sample(n=sample_size, shuffle=True).to_series().to_list()


def determine_species(df: pl.DataFrame, nr: int = 10, api_base_url: str = STRING_API_BASE_DEFAULT) -> int:
    """Determine species by sampling identifiers from a DataFrame and querying STRING API.

    Parameters:
        df: DataFrame with identifiers in the first column.
        nr: Number of identifiers to sample.
        api_base_url: STRING API base URL.

    Returns:
        The most frequent ncbiTaxonId found.

    Raises:
        ValueError: If the DataFrame is empty, has no valid identifiers,
                    or no taxon IDs could be retrieved.
    """
    sampled_ids = _sample_identifiers(df, nr)
    return determine_species_from_identifiers(sampled_ids, nr=len(sampled_ids), api_base_url=api_base_url)


def determine_species_from_identifiers(
    identifiers: list[str], nr: int = 10, api_base_url: str = STRING_API_BASE_DEFAULT
) -> int:
    """Determine species from a list of protein/gene identifiers.

    Parameters:
        identifiers: List of protein/gene identifiers.
        nr: Number of identifiers to sample (if list is larger).
        api_base_url: STRING API base URL.

    Returns:
        The most frequent ncbiTaxonId found.

    Raises:
        ValueError: If no taxon IDs could be retrieved.
    """
    import random

    sampled = random.sample(identifiers, min(nr, len(identifiers)))
    taxon_map = _fetch_ncbi_taxon_ids(sampled, api_base_url=api_base_url)
    taxon_ids = [tid for tid in taxon_map.values() if tid is not None]
    if not taxon_ids:
        raise ValueError("No taxon IDs returned from STRING API")
    return Counter(taxon_ids).most_common(1)[0][0]
