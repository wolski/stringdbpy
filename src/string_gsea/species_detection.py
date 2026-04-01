"""Species detection from FASTA OX fields in ZIP archives."""

import io
import re
import zipfile
from collections import Counter
from pathlib import Path

from loguru import logger

from string_gsea.gsea_config import STRING_API_BASE_DEFAULT
from string_gsea.models.gsea_models import RankListCollection
from string_gsea.string_api_client import determine_species_from_identifiers
from string_gsea.taxon_utils import TaxonUtils


def get_ox_fields(fasta_content: io.BytesIO) -> list[int]:
    """Extract OX= taxon IDs from FASTA header lines."""
    pattern = re.compile(r"OX=(\d+)")
    ox_values = []

    fasta_content = io.TextIOWrapper(fasta_content, encoding="utf-8")
    for line in fasta_content:
        if line.startswith(">"):
            match = pattern.search(line)
            if match:
                ox_values.append(int(match.group(1)))

    return ox_values


def get_species_from_oxes(zip_path: str) -> int:
    """Extract the most common species taxon ID from FASTA OX fields in a ZIP archive."""
    oxes = {}
    with zipfile.ZipFile(zip_path, "r") as z:
        fasta_files = [f for f in z.namelist() if f.endswith((".fas", ".fasta"))]
        logger.info(f"Found {len(fasta_files)} fasta files in {zip_path}")
        for file in fasta_files:
            logger.info(f"Processing file: {file}")
            with z.open(file) as f:
                file_bytes = f.read()
                ox = get_ox_fields(io.BytesIO(file_bytes))
                oxes[file] = ox

    merged = [item for sublist in oxes.values() for item in sublist]
    counter = Counter(merged)
    species = counter.most_common(1)[0][0]

    return int(species)


def get_species_taxon(
    zip_path: Path, rank_lists: RankListCollection, nr=10, api_base_url: str = STRING_API_BASE_DEFAULT
) -> int:
    """
    Determine the species taxon ID from the zip file.

    First tries to extract from FASTA OX fields, then validates against STRING.
    If not found in STRING, falls back to querying STRING API with protein identifiers.

    Args:
        zip_path: Path to the zip file containing rank/FASTA data
        rank_lists: Collection of RankList objects with protein identifiers
        nr: Number of identifiers to use for species determination (default: 10)
        api_base_url: STRING API base URL

    Returns:
        int: Valid STRING taxon ID
    """
    # Step 1: Extract raw taxon from FASTA OX fields
    raw_taxon = get_species_from_oxes(str(zip_path))

    # Step 2: Resolve to a STRING-DB supported species (walk taxonomy tree)
    string_taxon = TaxonUtils().get_organism_for_string(raw_taxon)
    if string_taxon is not None:
        logger.info(f"Species {string_taxon} resolved from FASTA OX fields")
        return string_taxon

    # Step 3: Fallback — query STRING API with sampled protein identifiers
    logger.info(f"OX-based taxon {raw_taxon} not in STRING, falling back to API lookup")
    first_rl = rank_lists.first()
    return determine_species_from_identifiers(list(first_rl.entries.keys()), nr=nr, api_base_url=api_base_url)
