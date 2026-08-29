"""Quarto rendering for completed STRING-GSEA result directories."""

import shutil
import subprocess
from pathlib import Path

from string_gsea.workflow.templates import TemplatePaths


def render_report(
    dataset_dir: Path,
    workunit_id: str,
    paths: TemplatePaths,
    done_file: Path,
) -> None:
    """Copy templates, render reports, clean temporary sources, and mark completion."""
    workunit = dataset_dir.resolve() / f"WU_{workunit_id}_GSEA"
    shutil.copy2(paths.vignettes / "GSEA_report.qmd", workunit / "GSEA_report.qmd")
    for name in ("index.qmd", "_fgcz-report.yml", "fgcz_header_quarto.html"):
        source = paths.templates / name
        if source.exists():
            shutil.copy2(source, workunit / name)
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
    for name in (
        "GSEA_report.qmd",
        "index.qmd",
        "_fgcz-report.yml",
        "fgcz_header_quarto.html",
        "index.rmarkdown",
    ):
        (workunit / name).unlink(missing_ok=True)
    for name in ("plots", "GSEA_report_files", "index_files"):
        directory = workunit / name
        if directory.exists():
            shutil.rmtree(directory)
    done_file.touch()
