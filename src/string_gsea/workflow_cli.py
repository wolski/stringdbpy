"""Composition roots for Snakemake, rendering, and packaging commands."""

import importlib.resources
import subprocess
import sys
from pathlib import Path

from cyclopts import App

from string_gsea.configuration import write_initial_configuration
from string_gsea.stringdb_adapters import RequestsApiKeyProvider
from string_gsea.workflow.packaging import package_results
from string_gsea.workflow.rendering import render_report
from string_gsea.workflow.templates import (
    OrderedTemplateLocator,
    RPackageTemplateLocator,
    TemplatePaths,
    WorkspaceTemplateLocator,
)

workflow_app = App(help="Run the STRING-GSEA Snakemake workflow.")
render_report_app = App(help="Render a STRING-GSEA report.")
render_only_app = App(help="Re-render an existing STRING-GSEA result directory.")
package_results_app = App(help="Package STRING-GSEA results for delivery.")


def _template_paths() -> TemplatePaths:
    repository_root = Path(__file__).resolve().parents[2]
    return OrderedTemplateLocator(
        (RPackageTemplateLocator(), WorkspaceTemplateLocator(repository_root))
    ).locate()


@workflow_app.default
def run_workflow(
    zip_path: str,
    workunit_id: str,
    out_dir: str = ".",
    *,
    which: str = "pep_2_no_imputed",
    fdr: float | None = None,
    cores: int = 1,
    dry_run: bool = False,
) -> None:
    """Run the packaged Snakemake workflow."""
    snakefile = importlib.resources.files("string_gsea.workflow") / "Snakefile"
    templates = _template_paths()
    config = [
        f"zip_path={zip_path}",
        f"workunit_id={workunit_id}",
        f"out_dir={out_dir}",
        f"which={which}",
        f"templates_dir={templates.templates}",
        f"vignettes_dir={templates.vignettes}",
    ]
    if fdr is not None:
        config.append(f"fdr={fdr}")
    command = [
        "snakemake",
        "-s",
        str(snakefile),
        "--cores",
        str(cores),
        "all",
        "--config",
        *config,
    ]
    if dry_run:
        command.append("-n")
    sys.exit(subprocess.run(command, check=False).returncode)


@workflow_app.command
def config() -> None:
    """Generate the STRING-GSEA configuration file."""
    write_initial_configuration(RequestsApiKeyProvider())


@render_report_app.default
def render_report_command(
    dataset_dir: str,
    workunit_id: str,
    vignettes_dir: str,
    templates_dir: str,
    done_file: str,
) -> None:
    """Render reports using explicitly supplied template directories."""
    render_report(
        Path(dataset_dir),
        workunit_id,
        TemplatePaths(templates=Path(templates_dir), vignettes=Path(vignettes_dir)),
        Path(done_file),
    )


@render_only_app.default
def render_only_command(dataset_dir: str, workunit_id: str) -> None:
    """Re-render reports using injected installed/workspace template locators."""
    result_directory = Path(dataset_dir) / f"WU_{workunit_id}_GSEA"
    render_report(
        Path(dataset_dir),
        workunit_id,
        _template_paths(),
        result_directory / "render_done.txt",
    )


@package_results_app.default
def package_results_command(output_base: str, workunit_id: str, outputs_yml_path: str) -> None:
    """Package results and write the B-Fabric manifest."""
    package_results(Path(output_base), workunit_id, Path(outputs_yml_path))
