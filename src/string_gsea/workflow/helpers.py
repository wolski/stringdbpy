"""Helper functions for the STRING-GSEA Snakemake workflow.

Functions are also exposed as CLI entry points so Snakemake shell: blocks
can invoke them without subprocess nesting issues from run: blocks.
"""

import logging
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

logger = logging.getLogger("string_gsea_workflow")


def _r_system_file(*args: str, package: str = "stringGSEAplot") -> Path | None:
    """Resolve installed R package paths via system.file()."""
    r_args = ", ".join(f"'{a}'" for a in args)
    r_code = f"cat(system.file({r_args}, package='{package}'))"
    try:
        result = subprocess.check_output(["Rscript", "-e", r_code], text=True, stderr=subprocess.DEVNULL).strip()
        if result:
            return Path(result)
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return None


def resolve_template_dirs() -> tuple[Path, Path]:
    """Resolve template and vignette paths.

    Returns (templates_dir, vignettes_dir):
    - templates_dir: index.qmd, _fgcz-report.yml, fgcz_header_quarto.html
    - vignettes_dir: GSEA_report.qmd
    """
    templates = _r_system_file("templates")
    vignettes = _r_system_file("doc")  # vignette sources installed here
    if templates and vignettes:
        return templates, vignettes

    # Local dev fallback
    r_pkg_dir = Path("../stringGSEAplot")
    return r_pkg_dir / "inst" / "templates", r_pkg_dir / "vignettes"


def render_report(
    dataset_dir: str,
    workunit_id: str,
    vignettes_dir: str,
    templates_dir: str,
    done_file: str,
) -> None:
    """Copy templates, render Quarto reports, clean up, and touch done file."""
    wu_dir = Path(dataset_dir).resolve() / f"WU_{workunit_id}_GSEA"
    json_file = f"WU{workunit_id}_gsea_result.json"

    # Copy QMD + support files into WU dir
    for name in ["GSEA_report.qmd"]:
        shutil.copy2(Path(vignettes_dir) / name, wu_dir / name)
    for name in ["index.qmd", "_fgcz-report.yml", "fgcz_header_quarto.html"]:
        src = Path(templates_dir) / name
        if src.exists():
            shutil.copy2(src, wu_dir / name)

    # Render GSEA report and index page
    subprocess.run(
        ["quarto", "render", "GSEA_report.qmd", "-P", f"json_path:{json_file}"],
        cwd=wu_dir,
        check=True,
    )
    subprocess.run(
        ["quarto", "render", "index.qmd", "-P", f"workunit_id:{workunit_id}", "-P", "package_dir:."],
        cwd=wu_dir,
        check=True,
    )

    # Clean up rendering artifacts
    for name in ["GSEA_report.qmd", "index.qmd", "_fgcz-report.yml", "fgcz_header_quarto.html", "index.rmarkdown"]:
        (wu_dir / name).unlink(missing_ok=True)
    for dirname in ["plots", "GSEA_report_files", "index_files"]:
        d = wu_dir / dirname
        if d.exists():
            shutil.rmtree(d)

    Path(done_file).touch()


def package_results(
    output_base: str,
    workunit_id: str,
    outputs_yml_path: str,
) -> None:
    """Zip the WU directory and write outputs.yml."""
    base = Path(output_base)
    wu_gsea = base / f"WU_{workunit_id}_GSEA"
    zip_path = base / f"WU_{workunit_id}_GSEA.zip"

    shutil.make_archive(
        str(zip_path).removesuffix(".zip"),
        "zip",
        root_dir=str(wu_gsea.parent),
        base_dir=wu_gsea.name,
    )
    logger.info("Package: %s", zip_path.resolve())

    outputs_data = {
        "outputs": [
            {
                "local_path": str(zip_path.resolve()),
                "store_entry_path": zip_path.name,
                "type": "bfabric_copy_resource",
            }
        ]
    }
    with open(outputs_yml_path, "w") as f:
        yaml.dump(outputs_data, f, default_flow_style=False)


def render_report_cli() -> None:
    """CLI entry point: render_report dataset_dir workunit_id vignettes_dir templates_dir done_file."""
    args = sys.argv[1:]
    if len(args) != 5:
        print("Usage: _string_gsea_render <dataset_dir> <workunit_id> <vignettes_dir> <templates_dir> <done_file>")
        sys.exit(1)
    render_report(*args)


def render_only_cli() -> None:
    """CLI entry point: re-render reports for an existing result dir without re-running STRING.

    Usage: _string_gsea_render_only <dataset_dir> <workunit_id>
    Resolves R package paths internally — useful inside Docker.
    """
    args = sys.argv[1:]
    if len(args) != 2:
        print("Usage: _string_gsea_render_only <dataset_dir> <workunit_id>")
        sys.exit(1)
    dataset_dir, workunit_id = args
    templates_dir, vignettes_dir = resolve_template_dirs()
    done_file = str(Path(dataset_dir) / f"WU_{workunit_id}_GSEA" / "render_done.txt")
    render_report(dataset_dir, workunit_id, str(vignettes_dir), str(templates_dir), done_file)


def package_results_cli() -> None:
    """CLI entry point: package_results output_base workunit_id outputs_yml_path."""
    args = sys.argv[1:]
    if len(args) != 3:
        print("Usage: _string_gsea_package <output_base> <workunit_id> <outputs_yml_path>")
        sys.exit(1)
    package_results(*args)
