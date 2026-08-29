# Changelog

## 0.0.9

- Standardize local and CI quality gates on the shared uv, Ruff, Pyright, deptry, pytest, and package-build workflow.
- Restore working source and wheel builds, including the packaged taxonomy mappings and Snakemake workflow.
- Replace the mutable GSEA builder/results flow and ORA monolith with typed `RunGSEA` and `RunORA` use cases using injected STRING-DB gateways.
- Preserve CLI commands and persisted GSEA/ORA artifacts while introducing polymorphic rank sources, analysis filters, species resolvers, and report-template locators.
- Enforce directed feature imports, strict Pyright for production and tests, Import Linter contracts, and 90% branch coverage in the shared quality checks.
- Publish and validate a packaged JSON Schema for ORA `enrichment_results.json` files.
