# Repository Guidelines

## Project Structure & Module Organization
- Root `*_cli.py` modules are composition roots. They construct adapters and inject them into the feature applications.
- `src/string_gsea/gsea/`, `ora/`, `taxonomy/`, and `workflow/` are mutually independent feature packages. Child packages never import root modules; domain models never import outward.
- `configuration.py` owns the TOML boundary and `stringdb_adapters.py` owns requests-based STRING-DB integrations. Consumer-owned protocols live in each feature's `ports.py`.
- Packaged resources live in `src/string_gsea/data/mappings/` and `src/string_gsea/workflow/Snakefile`.
- Tests reside in `tests/` with fixtures in `tests/data/`. Documentation and notebooks are under `docs/`.

## Build, Test, and Development Commands
- Synchronize the locked Python 3.13+ development environment with `make sync`.
- Run every merge-blocking gate with `make check`; use `make test`, `make lint`, `make typecheck`, `make deps`, or `make build` for narrower loops.
- Format source and tests with `make format`; verify formatting without edits with `make format-check`.
- Exercise CLI entry points with `--help`: `string_gsea_write_config`, `string_gsea_run`, `string_ora_run`, and `string_gsea_workflow`.

## Coding Style & Naming Conventions
- Format with Ruff (4-space indentation, 100-character lines). Type checking is strict; annotate public APIs and keep internal signatures honest.
- Modules and functions use `snake_case`; classes use `PascalCase`; constants are `UPPER_SNAKE`. Tests mirror target module names (for example, `test_session_yaml.py`).
- Depend on the smallest capability a function uses. Select variants through injected strategies, not type/mode discrimination or exception-driven fallback.
- Keep CLI commands, flags, output paths, JSON, YAML, and workflow resources compatible. Internal APIs may change without forwarding modules.

## Testing Guidelines
- Pytest is the primary framework; branch coverage must remain at or above the configured 90% gate.
- Name tests `test_*` and keep them deterministic—use fixtures in `tests/data/` rather than live STRING-DB calls. For integration-style checks (e.g., report rendering), prefer existing sample archives.
- Include targeted assertions for outputs (written files, TSV content, generated links) and durations for slow paths when practical.
- Production and test code must pass strict Pyright. Import Linter contracts and the Grimp architecture tests are merge-blocking.

## Commit & Pull Request Guidelines
- Use concise, imperative commit messages similar to the existing history (e.g., “Add CLAUDE.md…”, “Fix taxon problem”). Group related changes; avoid noisy churn.
- PRs should describe the change, note any new CLI flags or outputs, and list tests run (`make check` and any integration commands). Attach sample outputs or screenshots for report/notebook changes when feasible and link related issues.

## Configuration & Security Notes
- Sensitive values live in `$HOME/.config/string_gsea/config.toml` generated via `string_gsea_write_config`; never commit this file. Keep API keys and result ZIPs out of version control.
- Large result artifacts should be referenced, not stored. Prefer deterministic inputs in `tests/data/` for reproducible runs.
