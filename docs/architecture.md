# Directed architecture

STRING-GSEA keeps external behavior stable while making implementation dependencies explicit. Console commands and persisted files are compatibility boundaries; internal Python module paths are not.

## Dependency direction

```text
config_cli  gsea_cli  ora_cli  workflow_cli     composition roots
      \        |        |          /
       +--- stringdb_adapters ----+               requests implementations
                 |                                injected capabilities
      +----------+----------+----------+
      v          v          v          v
     gsea       ora      taxonomy   workflow       independent features
      |
 application -> ports / input / output
      |                    |
      +---------------> model                       inward domain values
```

The root CLI modules may combine feature packages. Feature packages do not import root modules or one another. Within GSEA, application and adapters point inward toward the model. `pyproject.toml` contains exhaustive Import Linter contracts, and `tests/test_architecture.py` adds graph and empty-initializer invariants.

## Use cases and injected boundaries

`RunGSEA.execute(RunGSEARequest) -> GSEAArtifacts` submits each typed rank list, waits for typed completed jobs, downloads optional external artifacts through `GSEAGateway`, parses results, and delegates persistence to `GSEAOutputWriter`.

`RunORA.execute(RunORARequest) -> ORAArtifacts` maps significant and background identifiers, requests enrichment through `ORAGateway`, and writes the established ORA directory through `ORAOutputWriter`.

The protocols are owned by their consumers. The root `RequestsStringDB` adapter implements the volatile HTTP capabilities. Polling receives monotonic-clock and sleep callables, so tests exercise success, failure, and timeout without waiting or making network calls. Requests exceptions are translated at the adapter boundary.

## Polymorphic selection

- `RankSource` implementations inspect an archive manifest and load XLSX or RNK data. The ordered source registry selects once; no exception controls fallback.
- `AnalysisPolicy` composes `RankFilter` objects for peptide-count and imputation rules.
- `OrderedSpeciesResolver` queries injected FASTA and identifier resolvers until one returns evidence. Exhaustion raises `SpeciesNotResolved`.
- `OrderedTemplateLocator` queries installed-package and workspace locators. Locators check executables, paths, and subprocess return codes explicitly.

Branches remain where they describe data or external state: malformed documents, empty identifiers, optional downloads, STRING status values, and missing files.

## Typed persistence boundaries

TOML, JSON, YAML, requests payloads, Polars data, and packaged mapping archives are validated before typed records are constructed. `SessionYaml` alone owns the established `outer~inner` encoding for tuple job keys. Characterization tests pin session YAML, GSEA JSON structure, output paths, archives, and `outputs.yml`.

## Verification

```bash
make check
uv run pre-commit run --all-files
```

The full check includes strict Pyright for `src/` and `tests/`, Ruff, Import Linter, deptry, deterministic pytest with at least 90% branch coverage, distribution builds, isolated wheel installation, packaged resource assertions, and console entry-point discovery. Network-dependent smoke and integration tests are marked and non-blocking.
