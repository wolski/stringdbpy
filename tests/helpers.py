"""Helper functions for the STRING-GSEA Snakemake workflow."""

import logging
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

import yaml

logger = logging.getLogger("string_gsea_workflow")


def get_datasets(datasets_dir: Path) -> list[str]:
    """Discover all datasets with params.yml that are not skipped."""
    datasets = []
    for params_file in sorted(datasets_dir.glob("*/params.yml")):
        with open(params_file) as f:
            params = yaml.safe_load(f)
        if not params.get("execution", {}).get("skip_run", False):
            datasets.append(params_file.parent.name)
    return datasets


@lru_cache
def load_params(dataset: str, datasets_dir: Path) -> dict:
    """Load params.yml for a dataset (cached)."""
    with open(datasets_dir / dataset / "params.yml") as f:
        return yaml.safe_load(f)


def get_ci_datasets(datasets: list[str], datasets_dir: Path) -> list[str]:
    """Filter datasets to those tagged with 'ci'."""
    return [
        d for d in datasets
        if "ci" in load_params(d, datasets_dir).get("metadata", {}).get("tags", [])
    ]


def _r_system_file(*args: str, package: str = "stringGSEAplot") -> Path | None:
    """Resolve installed R package paths via system.file()."""
    r_args = ", ".join(f'"{a}"' for a in args)
    cmd = f'Rscript -e "cat(system.file({r_args}, package=\\"{package}\\"))"'
    try:
        result = subprocess.check_output(
            cmd, shell=True, text=True, stderr=subprocess.DEVNULL
        ).strip()
        if result:
            return Path(result)
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return None


def resolve_template_dirs() -> tuple[Path, Path]:
    """Resolve template/vignette paths: installed R package first, local dev fallback."""
    installed = _r_system_file("templates")
    if installed:
        return installed, installed  # GSEA_report.qmd copied into templates at install

    r_pkg_dir = Path("../stringGSEAplot")
    return r_pkg_dir / "inst" / "templates", r_pkg_dir / "vignettes"


def run_gsea(
    zip_file: str,
    workunit_id: str,
    output_dir: str,
    from_rnk: bool,
    which: str | None,
    done_file: str,
) -> None:
    """Run string_gsea_run CLI and touch the done file."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    cmd = ["string_gsea_run", zip_file, workunit_id, output_dir]
    if from_rnk:
        cmd.append("--from-rnk")
    if which:
        cmd.extend(["--which", which])

    logger.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    Path(done_file).touch()


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
        shutil.copy2(Path(templates_dir) / name, wu_dir / name)

    # Render GSEA report and index page
    subprocess.run(
        ["quarto", "render", "GSEA_report.qmd", "-P", f"json_path:{json_file}"],
        cwd=wu_dir, check=True,
    )
    subprocess.run(
        ["quarto", "render", "index.qmd",
         "-P", f"workunit_id:{workunit_id}",
         "-P", "package_dir:."],
        cwd=wu_dir, check=True,
    )

    # Clean up rendering artifacts
    for name in ["GSEA_report.qmd", "index.qmd", "_fgcz-report.yml",
                  "fgcz_header_quarto.html", "index.rmarkdown"]:
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
        str(zip_path).removesuffix(".zip"), "zip",
        root_dir=str(wu_gsea.parent), base_dir=wu_gsea.name,
    )
    logger.info("Package: %s", zip_path.resolve())

    outputs_data = {"outputs": [{
        "local_path": str(zip_path.resolve()),
        "store_entry_path": zip_path.name,
        "type": "bfabric_copy_resource",
    }]}
    with open(outputs_yml_path, "w") as f:
        yaml.dump(outputs_data, f, default_flow_style=False)
