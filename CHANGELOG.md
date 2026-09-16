# Changelog

## 0.0.9

- Read ranks from the `AnnData.h5ad` a prolfquapp analysis writes, in preference to its `DE_*.xlsx` sheet. The artifact records which column is the contrast, the effect and the p-value, so a backend with its own column names — SAINTexpress reports `Bait`, `log2_EFCs`, `SaintScore`, `BFDR` — now ranks through the same code path as a linear model instead of failing on a missing `contrast` column.
- Rank scores from that artifact reproduce the `.rnk` files prolfquapp writes for the same analysis: signed `-log10(p)` where the model tests differences, and the effect size where it reports no p-value. This is asserted for every prolfqua modelling facade, SAINTexpress included.
- The `no_imputed` policies now drop estimates the model filled in rather than measured. On the AnnData path this reads the per-row estimate provenance; on the XLSX path the model-name pattern again matches current prolfquapp model names such as `lm_impute`, which it had silently stopped matching, leaving `pep_1_no_imputed` and `pep_2_no_imputed` doing nothing.
- The tabbed enrichment report and the result landing page now use the shared `fgczQuartoTemplate` look and feel — the FGCZ theme, banner and the top-right Find / Download / View-source toolbar — instead of a private copy of the report defaults. The report opens on an **Overview** tab with a visual abstract and an input summary, and closes on **Session Info** with report provenance and R session info.

- Standardize local and CI quality gates on the shared uv, Ruff, Pyright, deptry, pytest, and package-build workflow.
- Restore working source and wheel builds, including the packaged taxonomy mappings and Snakemake workflow.
- Replace the mutable GSEA builder/results flow and ORA monolith with typed `RunGSEA` and `RunORA` use cases using injected STRING-DB gateways.
- Preserve CLI commands and persisted GSEA/ORA artifacts while introducing polymorphic rank sources, analysis filters, species resolvers, and report-template locators.
- Enforce directed feature imports, strict Pyright for production and tests, Import Linter contracts, and 90% branch coverage in the shared quality checks.
- Publish and validate a packaged JSON Schema for ORA `enrichment_results.json` files.
