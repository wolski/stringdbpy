"""Run the STRING-GSEA Snakemake workflow."""

import importlib.resources
import subprocess
import sys

from cyclopts import App

app = App(help="Run the STRING-GSEA Snakemake workflow.")


@app.default
def main(
    zip_path: str,
    workunit_id: str,
    out_dir: str = ".",
    *,
    which: str = "pep_2_no_imputed",
    fdr: float | None = None,
    cores: int = 1,
    dry_run: bool = False,
) -> None:
    """Run STRING-GSEA analysis via Snakemake.

    Args:
        zip_path: Path to the input zip file.
        workunit_id: Identifier for this analysis run.
        out_dir: Path for output files.
        which: Analysis type (pep_1, pep_1_no_imputed, pep_2, pep_2_no_imputed, or none for RNK input).
        fdr: FDR threshold (overrides config.toml).
        cores: Number of cores for Snakemake.
        dry_run: Show what would be done without executing.
    """
    snakefile = importlib.resources.files("string_gsea.workflow") / "Snakefile"

    config_args = [
        f"zip_path={zip_path}",
        f"workunit_id={workunit_id}",
        f"out_dir={out_dir}",
        f"which={which}",
    ]
    if fdr is not None:
        config_args.append(f"fdr={fdr}")

    cmd = [
        "snakemake", "-s", str(snakefile),
        "--cores", str(cores),
        "all",
        "--config", *config_args,
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
