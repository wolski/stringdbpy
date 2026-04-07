"""Run the STRING-GSEA Snakemake workflow."""

import importlib.resources
import subprocess
import sys
from enum import Enum

from cyclopts import App

app = App(help="Run the STRING-GSEA Snakemake workflow.")


class Target(str, Enum):
    all = "all"  # GSEA + reports + packaging
    ci = "ci"  # CI-tagged datasets only
    clean = "clean"
    help = "help"


@app.default
def main(
    datasets_dir: str = "data/datasets",
    output_dir: str = "data/outputs",
    cores: int = 1,
    *,
    target: Target = Target.all,
    dry_run: bool = False,
) -> None:
    """Run STRING-GSEA batch processing via Snakemake.

    Args:
        datasets_dir: Path to directory containing dataset folders with params.yml.
        output_dir: Path for output files.
        cores: Number of cores for Snakemake.
        target: Workflow target (all, ci, clean, help).
        dry_run: Show what would be done without executing.
    """
    snakefile = importlib.resources.files("string_gsea.workflow") / "Snakefile"
    cmd = [
        "snakemake", "-s", str(snakefile),
        "--cores", str(cores),
        target.value,
        "--config", f"datasets_dir={datasets_dir}", f"output_dir={output_dir}",
    ]
    if dry_run:
        cmd.append("-n")

    result = subprocess.run(cmd)
    sys.exit(result.returncode)


@app.command
def config() -> None:
    """Generate STRING-GSEA configuration file."""
    from string_gsea.gsea_config import write_initial_configuration

    write_initial_configuration()
