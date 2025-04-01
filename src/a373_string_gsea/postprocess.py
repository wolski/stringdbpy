from pathlib import Path
import polars as pl
from pyexcelerate import Workbook
from typing import Dict

def write_xlsx(dfs: Dict[str, pl.DataFrame], filename: str) -> None:
    wb = Workbook()

    for sheet_name, df in dfs.items():
        print(f"Writing {sheet_name} with shape {df.shape}")

        # Build rows: header + data
        rows = [df.columns]

        for row in df.iter_rows(named=False):
            cleaned_row = [v if v is not None else None for v in row]
            rows.append(cleaned_row)

        wb.new_sheet(sheet_name, data=rows)

    wb.save(filename)

# Define the directory
directory = Path('../../WU_322935_GSEA')

# List all files (ignoring directories)
files = [f for f in directory.iterdir() if f.is_file()]
tsv_files = [f for f in files if f.suffix == '.tsv']

dfs = {f.name: pl.read_csv(f, separator = "\t") for f in tsv_files}

dfs_with_key = []
for key, df in dfs.items():
    # Add a new column "source" with the file name (key)
    df_with_key = df.with_columns(pl.lit(key).alias("contrast"))
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



# write comined_df to tsv
combined_df.write_csv("combined_df_longformat.tsv", separator="\t")
filtered_df = combined_df

columns = ["enrichmentScore", "genesInSet", "genesMapped", "direction","falseDiscoveryRate"]
pivot_dict = {}
for col in columns:
    pivot_dict[col] = filtered_df.pivot(
        values=col,
        index=["category", "termID", "termDescription","num_contrasts"],
        on="contrast"
    )

# merge the pivoted DataFrames
pivoted_dfs = []
for col, df in pivot_dict.items():
    # Prefix non-index columns with the variable name
    new_columns = [
        name if name in ["category", "termID", "termDescription", "num_contrasts"] else f"{col}_{name}"
        for name in df.columns
    ]
    df = df.rename(dict(zip(df.columns, new_columns)))
    pivoted_dfs.append(df)

# Join all pivoted DataFrames on the shared keys
merged_df = pivoted_dfs[0]
for df in pivoted_dfs[1:]:
    merged_df = merged_df.join(df, on=["category", "termID", "termDescription", "num_contrasts"], how="full")

x = pivot_dict["enrichmentScore"]
pivot_dict["DataUnfiltered"] = combined_df
pivot_dict['merged_df'] = merged_df

write_xlsx(pivot_dict, "string_gsea_results.xlsx")

print("done writing excel")

if __name__ == '__main__':
    species: int = 9606
    api_key = "b36F8oaRJwFZ"
    workunit_id = "322025"
