# Changelog

## 0.1.3

- The FGCZ-themed reports render again. The image shipped the report source without the template assets it includes, because the R package was installed with `build_vignettes=FALSE`, which bypasses the `vignettes/.install_extras` that ships them; the run produced every enrichment result and then died copying `_metadata.yml`.
- Rendering now goes through `stringGSEAplot::render_gsea_reports()`, which stages each report and the FGCZ assets via `fgczQuartoTemplate::fgcz_render()` exactly as prolfquapp does, so the theme, banner and toolbar come from the installed `fgczQuartoTemplate` at a pinned commit rather than from copies staged by this package.
- The image builds the R package from a tarball built with vignettes, which renders the packaged example report: a report that cannot render now fails the image build instead of a production run.

## 0.1.2

- Analyses from earlier prolfquapp versions rank again. prolfquapp has written `AnnData.h5ad` since 2.9.0 but only recorded the column roles from 2.10.0, and an archive holding a role-less artifact was claimed and then failed instead of being ranked from the `DE_*.xlsx` sheet beside it. An `.h5ad` now counts as a rank input only when its column roles can actually be read.
- The per-feature peptide count is read under either name prolfquapp has given it, `nrPeptides` or the earlier `nr_peptides`, so `pep_2` and `pep_2_no_imputed` work on results published before 2.10.4 as well. Requiring peptides of an analysis that records no count at all still fails, naming both.

## 0.1.1

- The minimum-peptides policies read `nrPeptides`, the one name prolfquapp now gives the per-feature peptide count. It used to be spelled `nr_peptides` in simulated analyses and `nrPeptides` in real ones, so `pep_2` and `pep_2_no_imputed` passed every fixture and then failed on a real B-Fabric analysis. Requires prolfquapp 2.10.4 or later.
- A `no_imputed` policy no longer fails on a SAINTexpress analysis. SAINTexpress imputes nothing, so its contrast table carried no `estimate_type` column at all and the filter had no provenance to read; prolfquasaint now stamps every row `observed`, and the policy is the no-op it should be. Every model is asserted to support all four policies.

## 0.1.0

- Make distribution smoke tests validate the artifacts built by the current version bump instead of looking for the preceding release's wheel.
- Read ranks from the `AnnData.h5ad` a prolfquapp analysis writes, in preference to its `DE_*.xlsx` sheet. The artifact records which column is the contrast, the effect and the p-value, so a backend with its own column names — SAINTexpress reports `Bait`, `log2_EFCs`, `SaintScore`, `BFDR` — now ranks through the same code path as a linear model instead of failing on a missing `contrast` column.
- Rank scores from that artifact reproduce the `.rnk` files prolfquapp writes for the same analysis: the recorded test statistic where the backend provides one, and the effect size where it reports no suitable score. This is asserted for every prolfqua modelling facade, SAINTexpress included.
- The `no_imputed` policies now drop estimates the model filled in rather than measured. On the AnnData path this reads the per-row estimate provenance and fails explicitly when that provenance is absent; on the XLSX path the model-name pattern again matches current prolfquapp model names such as `lm_impute`, which it had silently stopped matching, leaving `pep_1_no_imputed` and `pep_2_no_imputed` doing nothing.
- Selecting `--which none` now consumes the rank files already shipped in the archive instead of re-deriving an unfiltered ranking from its AnnData artifact.
- The tabbed enrichment report and the result landing page now use the shared `fgczQuartoTemplate` look and feel — the FGCZ theme, banner and the top-right Find / Download / View-source toolbar — instead of a private copy of the report defaults. The report opens on an **Overview** tab with a visual abstract and an input summary, and closes on **Session Info** with report provenance and R session info.

- Standardize local and CI quality gates on the shared uv, Ruff, Pyright, deptry, pytest, and package-build workflow; installing pre-commit now activates both the commit and push stages, with the 90% branch-coverage gate enforced before pushes.
- Restore working source and wheel builds, including the packaged taxonomy mappings and Snakemake workflow.
- Replace the mutable GSEA builder/results flow and ORA monolith with typed `RunGSEA` and `RunORA` use cases using injected STRING-DB gateways.
- Preserve CLI commands and persisted GSEA/ORA artifacts while introducing polymorphic rank sources, analysis filters, species resolvers, and report-template locators.
- Enforce directed feature imports, strict Pyright for production and tests, Import Linter contracts, and 90% branch coverage in the shared quality checks.
- Publish and validate a packaged JSON Schema for ORA `enrichment_results.json` files.
