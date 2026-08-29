# CLAUDE.md

## Project Overview

STRING-GSEA is a Python 3.13+ package and CLI for STRING-DB GSEA and ORA, species resolution, report rendering, and result packaging. The public compatibility boundary is the console commands and persisted artifacts; internal modules are intentionally architecture-driven.

## Development Commands

```bash
make sync          # synchronize the frozen uv development environment
make format        # Ruff format and safe fixes
make format-check  # verify formatting
make lint          # Ruff plus Import Linter contracts
make typecheck     # strict Pyright over src and tests
make deps          # deptry dependency validation
make test          # deterministic tests with >=90% branch coverage
make build         # sdist/wheel, metadata, resource, and entry-point checks
make check         # all merge-blocking checks
```

Network-dependent tests are separate: `make test-smoke` and `make test-integration`.

## Architecture

- `gsea_cli.py`, `ora_cli.py`, `workflow_cli.py`, and `config_cli.py` are composition roots.
- `stringdb_adapters.py` contains requests-backed implementations injected into consumer-owned ports.
- `gsea/`, `ora/`, `taxonomy/`, and `workflow/` are mutually independent. Children never import root modules.
- GSEA and ORA application classes express use cases (`RunGSEA`, `RunORA`). They receive narrow gateway capabilities in their constructors.
- Rank sources, rank filters, species resolvers, and template locators perform variant behavior polymorphically. Ordered registries select an implementation once.
- `configuration.py`, `gsea/session_yaml.py`, and the JSON model parsers validate untyped external data before constructing typed records.
- Package `__init__.py` files stay empty. Do not add forwarding modules or broad re-exports.

The executable rules are in `pyproject.toml` Import Linter contracts and `tests/test_architecture.py`. See `docs/architecture.md` for the dependency map.

## Compatibility Requirements

Preserve command names and arguments, `WU_{workunit_id}_GSEA` and `ORA_{workunit_id}` layouts, GSEA result JSON, session YAML (including `outer~inner` tuple-key encoding), `outputs.yml`, and packaged taxonomy/Snakefile resources.

The removed `StringGSEABuilder` and `StringGSEAResults` APIs are not compatibility surfaces. Do not restore them with shims; compose `RunGSEA` with a `GSEAGateway` instead.

## Running Commands

```bash
string_gsea_write_config --help
string_gsea_run --help
string_ora_run --help
string_gsea_workflow --help
```

Use `--which none` for RNK archives. The other supported analysis policies are `pep_1`, `pep_1_no_imputed`, `pep_2`, and `pep_2_no_imputed`.

The companion R package `stringGSEAplot` supplies report templates. Installed-package and workspace lookup are explicit injected locators; do not add exception-based fallback.
