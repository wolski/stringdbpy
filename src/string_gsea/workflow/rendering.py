"""Quarto rendering for completed STRING-GSEA result directories."""

import shutil
import subprocess
from pathlib import Path

from string_gsea.workflow.templates import TemplatePaths

# FGCZ Quarto template assets, vendored into stringGSEAplot beside the report
# sources by its `data-raw/sync_quarto_assets.R`. Quarto applies a file named
# `_metadata.yml` to every `.qmd` in its directory, so staging these next to the
# reports is what attaches the FGCZ theme, header and toolbar -- the reports
# themselves name no format.
FGCZ_ASSETS = (
    "_metadata.yml",
    "fgcz.scss",
    "fgcz_header_quarto.html",
    "fgcz-plot-finder.html",
)


def render_report(
    dataset_dir: Path,
    workunit_id: str,
    paths: TemplatePaths,
    done_file: Path,
) -> None:
    """Copy templates, render reports, clean temporary sources, and mark completion."""
    workunit = dataset_dir.resolve() / f"WU_{workunit_id}_GSEA"
    shutil.copy2(paths.vignettes / "GSEA_report.qmd", workunit / "GSEA_report.qmd")
    for name in FGCZ_ASSETS:
        shutil.copy2(paths.vignettes / name, workunit / name)
    shutil.copy2(paths.templates / "index.qmd", workunit / "index.qmd")
    subprocess.run(
        [
            "quarto",
            "render",
            "GSEA_report.qmd",
            "-P",
            f"json_path:WU{workunit_id}_gsea_result.json",
        ],
        cwd=workunit,
        check=True,
    )
    subprocess.run(
        [
            "quarto",
            "render",
            "index.qmd",
            "-P",
            f"workunit_id:{workunit_id}",
            "-P",
            "package_dir:.",
        ],
        cwd=workunit,
        check=True,
    )
    for name in ("GSEA_report.qmd", "index.qmd", "index.rmarkdown", *FGCZ_ASSETS):
        (workunit / name).unlink(missing_ok=True)
    for name in ("plots", "GSEA_report_files", "index_files"):
        directory = workunit / name
        if directory.exists():
            shutil.rmtree(directory)
    done_file.touch()
