"""Quarto rendering for completed STRING-GSEA result directories."""

import shutil
import subprocess
from pathlib import Path

from string_gsea.workflow.templates import TemplatePaths

# Rendering is delegated to the R package that owns the report sources, exactly
# as prolfquapp does: `protsea::render_gsea_reports()` stages each qmd
# and the FGCZ template assets through `fgczQuartoTemplate::fgcz_render()`, so
# the theme, banner and toolbar come from the installed `fgczQuartoTemplate`
# rather than from copies vendored into this package's own staging code.
# `Rscript -e` passes trailing words straight through to commandArgs(), so
# there is no `--args` separator here -- it would arrive as a fifth word.
_RENDER_CALL = (
    "args <- commandArgs(trailingOnly = TRUE); "
    "protsea::render_gsea_reports(args[[1]], args[[2]], args[[3]], args[[4]])"
)

# Quarto writes these beside the rendered HTML; the delivered archive carries
# the reports, not their intermediates.
RENDER_BYPRODUCTS = ("plots", "GSEA_report_files", "index_files")


def render_report(
    dataset_dir: Path,
    workunit_id: str,
    paths: TemplatePaths,
    done_file: Path,
) -> None:
    """Render the reports through protsea, then clean up and mark completion."""
    workunit = dataset_dir.resolve() / f"WU_{workunit_id}_GSEA"
    subprocess.run(
        [
            "Rscript",
            "-e",
            _RENDER_CALL,
            str(workunit),
            workunit_id,
            str(paths.vignettes),
            str(paths.templates),
        ],
        check=True,
    )
    for name in RENDER_BYPRODUCTS:
        directory = workunit / name
        if directory.exists():
            shutil.rmtree(directory)
    done_file.touch()
