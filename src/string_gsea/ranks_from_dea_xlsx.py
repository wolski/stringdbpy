import io
import zipfile
from pathlib import Path

import polars as pl
from loguru import logger

from string_gsea.models.gsea_models import RankList, RankListCollection


class DiffXLSX:
    def __init__(self, zip_path: Path):
        """
        Initialize the DiffXLSX class with a zip file path.

        Args:
            zip_path (Path): Path to the zip file containing the differential expression analysis Excel file.
        """
        self.zip_path = zip_path
        self.dea_df = DiffXLSX._read_diff_excel_from_zip(self.zip_path)

    @staticmethod
    def _read_diff_excel_from_zip(zip_path: Path) -> pl.DataFrame:
        """
        Read the differential expression Excel file from the zip archive.

        Args:
            zip_path (Path): Path to the zip file containing the differential expression analysis Excel file.

        Returns:
            pl.DataFrame: A Polars DataFrame containing the differential expression data.
        """
        # Open the zip archive
        with zipfile.ZipFile(zip_path, "r") as z:
            # List all files in the zip archive
            file_list = z.namelist()
            logger.debug(f"Files in zip: {file_list}")

            # Find the first .xlsx file in the archive
            xlsx_files = [file for file in file_list if file.endswith(".xlsx") and "DE_" in file]

            if not xlsx_files:
                raise ValueError("No xlsx file found in the zip archive")
            xlsx_filename = xlsx_files[0]
            logger.debug(f"Found xlsx file: {xlsx_filename}")

            # Open the xlsx file as a BytesIO stream
            with z.open(xlsx_filename) as f:
                xlsx_bytes = io.BytesIO(f.read())

                # Read all sheets from the Excel file into a dictionary
                df = pl.read_excel(xlsx_bytes, sheet_name="diff_exp_analysis")

        return df

    _VALID_ANALYSES = ("pep_1", "pep_1_no_imputed", "pep_2", "pep_2_no_imputed")

    @staticmethod
    def _get_ranks_by_contrast(df: pl.DataFrame) -> list[RankList]:
        """Build one RankList per contrast from a filtered DataFrame."""
        id_candidates = ["IDcolumn", "proteinname"]
        existing_id_cols = [col for col in id_candidates if col in df.columns]
        if not existing_id_cols:
            raise ValueError("No valid ID columns found")
        id_col = existing_id_cols[0]
        logger.info(f"Using ID column: {id_col}")

        result = []
        for contrast in df.select(pl.col("contrast")).unique().to_series().to_list():
            rank_df = df.filter(pl.col("contrast") == contrast).select(pl.col(id_col).alias("id"), pl.col("statistic"))
            result.append(RankList.from_polars(rank_df, contrast=str(contrast)))
        return result

    def rank_dict(
        self,
        which: str = "pep_2_no_imputed",
    ) -> RankListCollection:
        """
        Create rank lists for one analysis type from the differential expression data.

        Args:
            which: Analysis variant — one of "pep_1", "pep_1_no_imputed",
                   "pep_2", "pep_2_no_imputed".

        Returns:
            RankListCollection for the selected analysis type.
        """
        if which not in self._VALID_ANALYSES:
            raise ValueError(f"which must be one of: {', '.join(self._VALID_ANALYSES)}")

        filtered_df = self.dea_df
        if "no_imputed" in which:
            filtered_df = filtered_df.filter(~pl.col("modelName").str.contains("(?i)imputed"))
        if which.startswith("pep_2"):
            filtered_df = filtered_df.filter(pl.col("nrPeptides") > 1)

        return RankListCollection(
            analysis=which,
            rank_lists=self._get_ranks_by_contrast(filtered_df),
        )
