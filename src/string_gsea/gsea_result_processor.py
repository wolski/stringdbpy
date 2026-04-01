from pathlib import Path

import polars as pl
from loguru import logger
from pyexcelerate import Workbook

from string_gsea.models.gsea_models import GSEAResult


def _write_xlsx(dataframes: dict[str, pl.DataFrame], filename: Path) -> None:
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


def _to_wide(combined_df: pl.DataFrame, columns: list[str]) -> dict[str, pl.DataFrame]:
    pivot_dict = {}
    for col in columns:
        pivot_dict[col] = combined_df.pivot(
            values=col,
            index=["category", "termID", "termDescription", "num_contrasts"],
            on="contrast",
        )
    return pivot_dict


def _merge_pivoted_dfs(pivot_dict: dict[str, pl.DataFrame]) -> pl.DataFrame:
    join_cols = ["category", "termID", "termDescription", "num_contrasts"]
    pivoted_dfs = []
    for col, df in pivot_dict.items():
        # Prefix non-index columns with the variable name
        new_columns = [name if name in join_cols else f"{col}_{name}" for name in df.columns]
        df = df.rename(dict(zip(df.columns, new_columns, strict=True)))
        pivoted_dfs.append(df)

    merged_df = pivoted_dfs[0]
    for df in pivoted_dfs[1:]:
        merged_df = merged_df.join(df, on=join_cols, how="inner")
    return merged_df


def write_gsea_xlsx(gsea_result: GSEAResult, workunit_id: str, out_dir: Path) -> None:
    """Write long, pivoted, and merged XLSX reports from a GSEAResult."""
    combined_df = gsea_result.to_polars_long()

    # Long format
    _write_xlsx(
        {"longformat_df": combined_df},
        out_dir / f"WU{workunit_id}_string_gsea_results_long.xlsx",
    )

    # Pivoted
    pivot_columns = [
        "enrichmentScore",
        "genesInSet",
        "genesMapped",
        "directionNR",
        "falseDiscoveryRate",
    ]
    pivoted_dict = _to_wide(combined_df, pivot_columns)
    _write_xlsx(
        pivoted_dict,
        out_dir / f"WU{workunit_id}_string_gsea_results_pivoted.xlsx",
    )

    # Merged
    merged_df = _merge_pivoted_dfs(pivoted_dict)
    _write_xlsx(
        {"merged_df": merged_df},
        out_dir / f"WU{workunit_id}_string_gsea_results_merged.xlsx",
    )
