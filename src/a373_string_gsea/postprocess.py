from pathlib import Path
import polars as pl
from pyexcelerate import Workbook
from typing import Dict, List
from loguru import logger


def write_xlsx(dataframes: Dict[str, pl.DataFrame], filename: Path) -> None:
    wb = Workbook()
    for sheet_name, dataframe in dataframes.items():
        logger.info(f"Writing {sheet_name} with shape {dataframe.shape}")
        rows = [dataframe.columns]
        for row in dataframe.iter_rows(named=False):
            cleaned_row = list(row)
            rows.append(cleaned_row)

        wb.new_sheet(sheet_name, data=rows)
    wb.save(filename)

    logger.info(f"Wrote {filename}")


def results_to_dataframe(tsv_files: List[Path]) -> pl.DataFrame:
    dfs_with_key = []
    dfs = {f.name: pl.read_csv(f, separator="\t") for f in tsv_files}
    for key, dataframe in dfs.items():
        # Add a new column "source" with the file name (key)
        df_with_key = dataframe.with_columns(pl.lit(key).alias("contrast"))
        # Reorder columns so that "source" is first
        cols = df_with_key.columns
        new_order = ['contrast'] + [col for col in cols if col != 'contrast']
        df_with_key = df_with_key.select(new_order)
        dfs_with_key.append(df_with_key)
    # Concatenate all DataFrames into a single DataFrame
    combined_df = pl.concat(dfs_with_key)
    combined_df = combined_df.with_columns(
        pl.when(pl.col("direction") == "top").then(1)
        .when(pl.col("direction") == "bottom").then(-1)
        .otherwise(0).alias("directionNR")
    )
    grouped = combined_df.group_by(["category", "termID"]).agg(pl.count("contrast").alias("num_contrasts"))
    combined_df = combined_df.join(grouped, on=["category", "termID"], how="inner")
    return combined_df


def to_wide(combined_df: pl.DataFrame, columns: List[str]) -> Dict[str, pl.DataFrame]:
    pivot_dict = {}
    for col in columns:
        pivot_dict[col] = combined_df.pivot(
            values=col,
            index=["category", "termID", "termDescription", "num_contrasts"],
            on="contrast"
        )
    return pivot_dict


def merge_pivoted_dfs(pivot_dict: Dict[str, pl.DataFrame]) -> pl.DataFrame:
    join_cols = ["category", "termID", "termDescription", "num_contrasts"]
    pivoted_dfs = []
    for col, df in pivot_dict.items():
        # Prefix non-index columns with the variable name
        new_columns = [
            name if name in join_cols else f"{col}_{name}"
            for name in df.columns
        ]
        df = df.rename(dict(zip(df.columns, new_columns)))
        pivoted_dfs.append(df)

    merged_df = pivoted_dfs[0]
    for df in pivoted_dfs[1:]:
        merged_df = merged_df.join(df, on=join_cols, how="inner")
    return merged_df

def list_subfolders(path: Path):
    return [p for p in path.iterdir() if p.is_dir()]

def _result_to_xlsx(directory: Path, workunit_id: str) -> None:
    tsv_files = list(directory.glob("*.tsv"))

    logger.info(f"Found {len(tsv_files)} tsv files")
    if len(tsv_files) == 0:
        logger.error("No tsv files found")
        return

    combined_df = results_to_dataframe(tsv_files)
    pivoted_dict = to_wide(combined_df,
                           ["enrichmentScore", "genesInSet", "genesMapped", "directionNR", "falseDiscoveryRate"])

    merged_df = merge_pivoted_dfs(pivoted_dict)
    write_xlsx(pivoted_dict, directory / f"WU{workunit_id}_string_gsea_results_pivoted.xlsx")
    logger.info("written pivoted XLSX")
    merged_df = {"merged_df": merged_df}
    write_xlsx(merged_df, directory / f"WU{workunit_id}_string_gsea_results_merged.xlsx")
    logger.info("written merged XLSX")

    longformat_df = {"longformat_df": combined_df}
    write_xlsx(longformat_df, directory / f"WU{workunit_id}_string_gsea_results_long.xlsx")
    logger.info("written merged XLSX")

def result_to_xlsx(w_directory:Path, workunit_id: str) -> None:
    directories = list_subfolders(w_directory)
    for directory in directories:
        _result_to_xlsx(directory, workunit_id)

if __name__ == '__main__':
    # Define the directory
    w_directory = Path('./WU_234567_GSEA')
    result_to_xlsx(w_directory, "234567")
    print("done writing excel")
