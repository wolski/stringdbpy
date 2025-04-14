import zipfile
import io
from pathlib import Path

import polars as pl
from loguru import logger

from a373_string_gsea.run_string_gsea import extract_workunit_id_from_file, find_zip_files, get_species_from_oxes
from a373_string_gsea.stringgsea import StringGSEA


def read_diff_excel_from_zip(zip_path: str) -> pl.DataFrame:
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


def get_ranks_by_contrast(df: pl.DataFrame, id_col: str = "IDcolumn", rank_col: str = "statistic",
                          prefix: str = "") -> dict:
    df_dict = {
        prefix + str(contrast): df.filter(pl.col("contrast") == contrast).select([id_col, rank_col])
        for contrast in df.select(pl.col("contrast")).unique().to_series().to_list()
    }
    return df_dict


def make_rank_lists(dea_df: pl.DataFrame, id_col: str = "IDcolumn", rank_col: str = "statistic") -> dict:
    contrast_pep_1 = get_ranks_by_contrast(dea_df, id_col, rank_col)
    df_pep1_no_imp = dea_df.filter(~pl.col("modelName").str.contains("(?i)imputed"))
    contrast_pep1_no_imp = get_ranks_by_contrast(df_pep1_no_imp, id_col, rank_col)
    df_pep2 = dea_df.filter(pl.col("nrPeptides") > 1)
    contrast_pep2 = get_ranks_by_contrast(df_pep2, id_col, rank_col)
    df_pep2_no_imp = df_pep2.filter(~pl.col("modelName").str.contains("(?i)imputed"))
    contrast_pep2_no_imp = get_ranks_by_contrast(df_pep2_no_imp, id_col, rank_col)
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
    species: int
    api_key = "b36F8oaRJwFZ"
    workunit_id: str = "1234567"

    if Path("params.yml").exists():
        workunit_id = str(extract_workunit_id_from_file("params.yml"))
    else:
        workunit_id = 234567

    logger.info(f"Workunit ID: {workunit_id}")
    fdr: float = 0.25
    zip_path = "../../testing_examples/fromDEA/2790797.zip"
    species = get_species_from_oxes(zip_path)

    df = read_diff_excel_from_zip(zip_path)
    rank_files = make_rank_lists(df)

    gsea = StringGSEA(api_key, workunit_id, rank_files, species, fdr)
    gsea.submit()
    logger.info(f"Job submitted successfully.{gsea.res_job_id}")
    gsea.pull_results()
    logger.info("got results")
    gsea.serialize_results()
    gsea.write_links()
    gsea.write_gsea_results()
    result_dir = gsea.write_rank_files()
    gsea.zip_folder(result_dir)

