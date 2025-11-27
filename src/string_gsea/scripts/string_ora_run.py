from cyclopts import App
from pathlib import Path
from loguru import logger
import requests
import re
import json
import polars as pl
from collections import Counter

from string_gsea.gsea_config import get_configuration
from string_gsea.get_species import TaxonUtils

app = App()

STRING_API_BASE = "https://version-12-0.string-db.org/api"


def get_species_from_fasta(fasta_path: Path) -> int:
    """
    Extract species taxon ID from FASTA file OX= fields.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        int: NCBI Taxon ID validated against STRING-DB
    """
    pattern = re.compile(r'OX=(\d+)')
    ox_values = []

    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                match = pattern.search(line)
                if match:
                    ox_values.append(int(match.group(1)))

    if not ox_values:
        raise ValueError(f"No OX= fields found in FASTA file: {fasta_path}")

    counter = Counter(ox_values)
    taxon = counter.most_common(1)[0][0]

    # Validate against STRING-DB supported species
    taxon_utils = TaxonUtils()
    valid_taxon = taxon_utils.get_organism_for_string(taxon)

    if valid_taxon is None:
        raise ValueError(f"Taxon {taxon} is not supported by STRING-DB")

    logger.info(f"Species taxon ID: {valid_taxon}")
    return valid_taxon


def read_id_list(filepath: Path) -> list[str]:
    """
    Read protein/gene identifiers from a file (one per line).

    Args:
        filepath: Path to text file with one identifier per line

    Returns:
        list[str]: List of identifiers
    """
    with open(filepath, 'r') as f:
        ids = [line.strip() for line in f if line.strip()]
    logger.info(f"Read {len(ids)} identifiers from {filepath}")
    return ids


def map_to_string_ids(
    identifiers: list[str],
    species: int,
    caller_identity: str,
    batch_size: int = 500
) -> dict[str, str]:
    """
    Map identifiers to STRING IDs using the STRING API.

    Args:
        identifiers: List of protein/gene identifiers
        species: NCBI Taxon ID
        caller_identity: Caller identity string
        batch_size: Number of IDs per batch request

    Returns:
        dict: Mapping from input identifier to STRING ID
    """
    url = f"{STRING_API_BASE}/json/get_string_ids"
    mapping = {}

    # Process in batches to avoid URL length limits
    for i in range(0, len(identifiers), batch_size):
        batch = identifiers[i:i + batch_size]
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "caller_identity": caller_identity,
        }

        logger.info(f"Mapping batch {i // batch_size + 1}: {len(batch)} identifiers")
        resp = requests.get(url, params=params)
        resp.raise_for_status()

        data = resp.json()
        for entry in data:
            query_item = entry.get("queryItem")
            string_id = entry.get("stringId")
            if query_item and string_id:
                mapping[query_item] = string_id

    logger.info(f"Mapped {len(mapping)} of {len(identifiers)} identifiers to STRING IDs")
    return mapping


def get_string_link(
    identifiers: list[str],
    species: int,
    caller_identity: str
) -> str:
    """
    Get a link to view the network on STRING-DB web interface.

    Args:
        identifiers: List of STRING IDs
        species: NCBI Taxon ID
        caller_identity: Caller identity string

    Returns:
        str: URL to STRING-DB network visualization
    """
    url = f"{STRING_API_BASE}/json/get_link"
    params = {
        "identifiers": "\r".join(identifiers),
        "species": species,
        "caller_identity": caller_identity,
    }

    # Use POST to avoid URL length limits with many identifiers
    resp = requests.post(url, data=params)
    resp.raise_for_status()

    data = resp.json()
    if isinstance(data, list) and len(data) > 0:
        return data[0]
    return None


def run_enrichment(
    identifiers: list[str],
    background_string_ids: list[str],
    species: int,
    caller_identity: str
) -> list[dict]:
    """
    Run functional enrichment analysis via STRING-DB API with custom background.

    Args:
        identifiers: List of STRING IDs for significant genes
        background_string_ids: List of STRING IDs for background
        species: NCBI Taxon ID
        caller_identity: Caller identity string

    Returns:
        list[dict]: Enrichment results
    """
    url = f"{STRING_API_BASE}/json/enrichment"
    params = {
        "identifiers": "\r".join(identifiers),
        "background_string_identifiers": "\r".join(background_string_ids),
        "species": species,
        "caller_identity": caller_identity,
    }

    logger.info(f"Running enrichment for {len(identifiers)} identifiers with {len(background_string_ids)} background")
    resp = requests.post(url, data=params)
    resp.raise_for_status()

    data = resp.json()

    # Check for error response
    if isinstance(data, list) and len(data) > 0:
        if data[0].get('error'):
            msg = data[0].get('message', 'Unknown error')
            raise RuntimeError(f"STRING-DB error: {msg}")

    logger.info(f"Received {len(data)} enrichment terms")
    return data


def save_results(results: list[dict], out_dir: Path) -> dict:
    """
    Save enrichment results to JSON and TSV files.

    Args:
        results: Enrichment results from API
        out_dir: Output directory

    Returns:
        dict: Paths to saved files
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}

    # Save JSON
    json_path = out_dir / "enrichment_results.json"
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    paths['json'] = json_path
    logger.info(f"Saved JSON: {json_path}")

    # Convert to TSV using polars
    if results:
        df = pl.DataFrame(results)
        # Convert list columns to comma-separated strings for TSV compatibility
        for col in df.columns:
            if df[col].dtype == pl.List:
                df = df.with_columns(
                    pl.col(col).list.join(",").alias(col)
                )
        tsv_path = out_dir / "enrichment_results.tsv"
        df.write_csv(tsv_path, separator='\t')
        paths['tsv'] = tsv_path
        logger.info(f"Saved TSV: {tsv_path}")

    return paths


@app.default()
def string_ora_run(
    significant: str,
    background: str,
    fasta: str,
    out_dir: str = ".",
    workunit_id: str = "ORA"
):
    """
    Run STRING ORA (Over-Representation Analysis) on the provided gene lists.

    Args:
        significant: Path to significant.txt (one protein/gene ID per line)
        background: Path to background.txt (one protein/gene ID per line)
        fasta: Path to FASTA file (used to determine organism via OX= field)
        out_dir: Base directory for output files. Defaults to current directory.
        workunit_id: Identifier for this analysis run. Defaults to "ORA".
    """
    # Convert paths
    significant_path = Path(significant)
    background_path = Path(background)
    fasta_path = Path(fasta)
    base_dir = Path(out_dir)

    # Validate inputs
    if not significant_path.exists():
        raise FileNotFoundError(f"Significant file not found: {significant_path}")
    if not background_path.exists():
        raise FileNotFoundError(f"Background file not found: {background_path}")
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    # Load configuration
    config = get_configuration()

    # Determine species from FASTA
    species = get_species_from_fasta(fasta_path)

    # Read gene lists
    sig_ids = read_id_list(significant_path)
    bg_ids = read_id_list(background_path)

    if not sig_ids:
        raise ValueError("No identifiers in significant file")
    if not bg_ids:
        raise ValueError("No identifiers in background file")

    # Create output directory
    result_dir = base_dir / f"ORA_{workunit_id}"
    result_dir.mkdir(parents=True, exist_ok=True)

    # Map identifiers to STRING IDs
    logger.info("Mapping significant identifiers to STRING IDs...")
    sig_mapping = map_to_string_ids(sig_ids, species, config.caller_identity)

    logger.info("Mapping background identifiers to STRING IDs...")
    bg_mapping = map_to_string_ids(bg_ids, species, config.caller_identity)

    # Get mapped STRING IDs
    sig_string_ids = list(sig_mapping.values())
    bg_string_ids = list(bg_mapping.values())

    if not sig_string_ids:
        raise ValueError("No significant identifiers could be mapped to STRING IDs")
    if not bg_string_ids:
        raise ValueError("No background identifiers could be mapped to STRING IDs")

    logger.info(f"Mapped {len(sig_string_ids)} significant and {len(bg_string_ids)} background to STRING IDs")

    # Run enrichment with background
    results = run_enrichment(
        identifiers=sig_string_ids,
        background_string_ids=bg_string_ids,
        species=species,
        caller_identity=config.caller_identity
    )

    # Save results
    paths = save_results(results, result_dir)

    # Get and save STRING-DB link
    logger.info("Getting STRING-DB network link...")
    string_link = get_string_link(
        identifiers=sig_string_ids,
        species=species,
        caller_identity=config.caller_identity
    )
    if string_link:
        links_path = result_dir / "links.txt"
        with open(links_path, 'w') as f:
            f.write(f"STRING-DB Network: {string_link}\n")
        paths['links'] = links_path
        logger.info(f"Saved link: {links_path}")

    # Save mapping info
    mapping_path = result_dir / "id_mapping.json"
    with open(mapping_path, 'w') as f:
        json.dump({
            "significant": sig_mapping,
            "background": bg_mapping,
            "unmapped_significant": [x for x in sig_ids if x not in sig_mapping],
            "unmapped_background": [x for x in bg_ids if x not in bg_mapping],
        }, f, indent=2)
    paths['mapping'] = mapping_path
    logger.info(f"Saved ID mapping: {mapping_path}")

    logger.info(f"ORA analysis complete. Results in: {result_dir}")
    return paths


if __name__ == '__main__':
    app()
