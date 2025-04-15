import zipfile
import io
from pathlib import Path

import polars as pl
from loguru import logger




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
        with zipfile.ZipFile(zip_path, 'r') as z:
            # List all files in the zip archive
            file_list = z.namelist()
            print("Files in zip:", file_list)

            # Find the first .xlsx file in the archive
            xlsx_files = [file for file in file_list if file.endswith('.xlsx') and 'DE_' in file]

            if not xlsx_files:
                raise ValueError("No xlsx file found in the zip archive")
            xlsx_filename = xlsx_files[0]
            print("Found xlsx file:", xlsx_filename)

            # Open the xlsx file as a BytesIO stream
            with z.open(xlsx_filename) as f:
                xlsx_bytes = io.BytesIO(f.read())

                # Read all sheets from the Excel file into a dictionary
                df = pl.read_excel(xlsx_bytes, sheet_name="diff_exp_analysis")

        return df


    @staticmethod
    def _get_ranks_by_contrast(df: pl.DataFrame, id_col: str = "IDcolumn", rank_col: str = "statistic",
                              prefix: str = "") -> dict:
        """
        Get rank dataframes by contrast.
        
        Args:
            df (pl.DataFrame): The input DataFrame.
            id_col (str): The column name for identifiers.
            rank_col (str): The column name for rank values.
            prefix (str): A prefix to add to the contrast names.
            
        Returns:
            dict: A dictionary mapping contrast names to rank DataFrames.
        """
        df_dict = {
            prefix + str(contrast): df.filter(pl.col("contrast") == contrast).select([id_col, rank_col])
            for contrast in df.select(pl.col("contrast")).unique().to_series().to_list()
        }
        return df_dict


    def rank_dict(self, id_col: str = "IDcolumn", rank_col: str = "statistic") -> dict:
        """
        Create rank lists from the differential expression DataFrame.
        
        Args:
            id_col (str): The column name for identifiers.
            rank_col (str): The column name for rank values.
            
        Returns:
            dict: A dictionary mapping (peptide_type, contrast) tuples to rank DataFrames.
        """
        contrast_pep_1 = self._get_ranks_by_contrast(self.dea_df, id_col, rank_col)
        df_pep1_no_imp = self.dea_df.filter(~pl.col("modelName").str.contains("(?i)imputed"))
        contrast_pep1_no_imp = self._get_ranks_by_contrast(df_pep1_no_imp, id_col, rank_col)
        df_pep2 = self.dea_df.filter(pl.col("nrPeptides") > 1)
        contrast_pep2 = self._get_ranks_by_contrast(df_pep2, id_col, rank_col)
        df_pep2_no_imp = df_pep2.filter(~pl.col("modelName").str.contains("(?i)imputed"))
        contrast_pep2_no_imp = self._get_ranks_by_contrast(df_pep2_no_imp, id_col, rank_col)
        combined_dict = {
            "pep_1": contrast_pep_1,
            "pep_1_no_imputed": contrast_pep1_no_imp,
            "pep_2": contrast_pep2,
            "pep_2_no_imputed": contrast_pep2_no_imp}
        flattened = {}
        for outer_key, subdict in combined_dict.items():
            for inner_key, rank_df in subdict.items():
                flattened[(outer_key, inner_key)] = rank_df

        return flattened


if __name__ == '__main__':

    zip_path = Path(__file__).parent.parent.parent / "tests/data/DE_mouse_fasta_xlsx.zip"
    logger.info(f"Zip path: {zip_path}")
    df_xlsx = DiffXLSX(zip_path)
    rank_files = df_xlsx.rank_dict()


