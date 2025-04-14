import re
import zipfile
import io
import polars as pl
from pathlib import Path
from loguru import logger


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
                filename_only = Path(file).stem
                dataframes[("from_rnk", filename_only)] = df
    return dataframes



def find_zip_files(path: Path = Path(".")):
    pattern = re.compile(r'^(\d{7}.*|DEA_[^/]+)\.zip$')

    candidates = list(path.glob("*.zip"))
    if not candidates:
        logger.error("No zip files found")
        raise FileNotFoundError(f"No zip files found in {path.resolve()}")

    # Filter files whose names match the pattern
    files = [f for f in candidates if pattern.match(f.name)]
    if not files:
        logger.error("No zip files match the pattern.")
        raise FileNotFoundError("No zip files found matching the pattern.")

    return files
