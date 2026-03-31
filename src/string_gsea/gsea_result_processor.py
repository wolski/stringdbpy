import tempfile
from pathlib import Path

import polars as pl
from loguru import logger
from pyexcelerate import Workbook

from string_gsea.models.gsea_models import parse_gsea_tsv_dir


class GSEAResultProcessor:
    @staticmethod
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

    @staticmethod
    def _to_wide(combined_df: pl.DataFrame, columns: list[str]) -> dict[str, pl.DataFrame]:
        pivot_dict = {}
        for col in columns:
            pivot_dict[col] = combined_df.pivot(
                values=col,
                index=["category", "termID", "termDescription", "num_contrasts"],
                on="contrast",
            )
        return pivot_dict

    @staticmethod
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

    @staticmethod
    def _list_subfolders(path: Path):
        return [p for p in path.iterdir() if p.is_dir()]

    @staticmethod
    def _process_directory(directory: Path, workunit_id: str) -> None:
        try:
            gsea_result = parse_gsea_tsv_dir(directory)
        except FileNotFoundError:
            logger.error("No tsv files found")
            return

        # Write JSON
        json_path = directory / f"WU{workunit_id}_gsea_result.json"
        gsea_result.to_json(json_path)
        logger.info(f"Written JSON: {json_path}")

        # Write XLSX (long, pivoted, merged)
        combined_df = gsea_result.to_polars_long()
        pivoted_dict = GSEAResultProcessor._to_wide(
            combined_df,
            [
                "enrichmentScore",
                "genesInSet",
                "genesMapped",
                "directionNR",
                "falseDiscoveryRate",
            ],
        )

        merged_df = GSEAResultProcessor._merge_pivoted_dfs(pivoted_dict)
        GSEAResultProcessor._write_xlsx(
            pivoted_dict,
            directory / f"WU{workunit_id}_string_gsea_results_pivoted.xlsx",
        )
        logger.info("Written pivoted XLSX")
        merged_df = {"merged_df": merged_df}
        GSEAResultProcessor._write_xlsx(merged_df, directory / f"WU{workunit_id}_string_gsea_results_merged.xlsx")
        logger.info("Written merged XLSX")

        longformat_df = {"longformat_df": combined_df}
        GSEAResultProcessor._write_xlsx(longformat_df, directory / f"WU{workunit_id}_string_gsea_results_long.xlsx")
        logger.info("Written long XLSX")

    @classmethod
    def process_results(cls, work_directory: Path, workunit_id: str) -> None:
        directories = cls._list_subfolders(work_directory)
        for directory in directories:
            cls._process_directory(directory, workunit_id)


if __name__ == "__main__":
    # Define the directory
    current_dir = Path(__file__).parent.parent.parent
    w_directory = current_dir / "tests" / "data" / "dummy_out" / "WU_abcd_GSEA"
    # check if the directory exists
    if not w_directory.exists():
        raise FileNotFoundError(f"Directory {w_directory} does not exist")
    # get proper temp directory
    GSEAResultProcessor.process_results(w_directory, "234567")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        print("done writing excel")
