# Repository Guidelines

## Project Structure & Module Organization
- `src/string_gsea/` holds the Python package. Key modules include the builder (`string_gsea_builder.py`), session management (`gsea_session.py`), results handling (`string_gsea_results.py`), configuration (`gsea_config.py`/`config.py`), plotting/network utilities, and CLI entry points under `scripts/`.
- Assets and templates live in `src/string_gsea/data/` (mapping ZIPs) and `src/string_gsea/docs/` (Quarto report sources). External helper scripts are under `scripts/`.
- Tests reside in `tests/` with fixtures in `tests/data/`. Documentation and notebooks are under `docs/` (see `docs/marimo_notebooks/` for Marimo examples).

## Build, Test, and Development Commands
- Install for development with test extras: `uv pip install -e ".[test]"` (Python 3.13+). Use `uv pip install -e ".[notebooks]"` for Marimo work.
- Run the full suite via nox (uses UV backend): `nox -s test`. Direct pytest is fine for local loops: `pytest --durations=50 tests -k <pattern>`.
- Exercise CLI entry points: `string_gsea_write_config --help`, `string_gsea_run --help`, and `string_gsea_render_reports --help`. Internal script smoke tests: `nox -s run-internal-scripts` (or pass module paths to narrow scope).

## Coding Style & Naming Conventions
- Format with Black defaults (4-space indentation, 88-char lines). Prefer type hints and lightweight docstrings for public APIs.
- Modules and functions use `snake_case`; classes use `PascalCase`; constants are `UPPER_SNAKE`. Tests mirror target module names (e.g., `test_gsea_session.py`).
- Keep CLI surfaces consistent with existing flags and parameter names; avoid breaking changes without deprecation notes.

## Testing Guidelines
- Pytest is the primary framework; coverage is configured in `pyproject.toml`. Aim to maintain or raise coverage when modifying core workflows.
- Name tests `test_*` and keep them deterministic—use fixtures in `tests/data/` rather than live STRING-DB calls. For integration-style checks (e.g., report rendering), prefer existing sample archives.
- Include targeted assertions for outputs (written files, TSV content, generated links) and durations for slow paths when practical.

## Commit & Pull Request Guidelines
- Use concise, imperative commit messages similar to the existing history (e.g., “Add CLAUDE.md…”, “Fix taxon problem”). Group related changes; avoid noisy churn.
- PRs should describe the change, note any new CLI flags or outputs, and list tests run (`nox -s test`, targeted `pytest` commands). Attach sample outputs or screenshots for report/notebook changes when feasible and link related issues.

## Configuration & Security Notes
- Sensitive values live in `$HOME/.config/string_gsea/config.toml` generated via `string_gsea_write_config`; never commit this file. Keep API keys and result ZIPs out of version control.
- Large result artifacts should be referenced, not stored. Prefer deterministic inputs in `tests/data/` for reproducible runs.
