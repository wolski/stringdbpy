"""CLI script for rendering Marimo-based GSEA reports."""

import importlib.resources
import shlex
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from cyclopts import App
from loguru import logger

from string_gsea.string_gsea_results import StringGSEAResults

app = App(
    help="Render Marimo reports as static HTML files."
)


def execute_marimo_export(
    notebook_path: Path,
    output_file: Path,
    data_file: Path,
    links_file: Path,
    fdr_threshold: float = 0.05,
    genes_mapped_threshold: int = 10,
    max_terms: int = 100,
) -> bool:
    """
    Export a Marimo notebook to static HTML.

    Args:
        notebook_path: Path to the Marimo notebook (.py file)
        output_file: Path for the output HTML file
        data_file: Path to the XLSX results file
        links_file: Path to the links file
        fdr_threshold: FDR threshold for filtering
        genes_mapped_threshold: Minimum genes mapped threshold
        max_terms: Maximum terms per category/contrast to render

    Returns:
        True if export succeeded, False otherwise
    """
    # Use sys.executable to find marimo in the same environment
    marimo_path = Path(sys.executable).parent / "marimo"
    cmd = [
        str(marimo_path),
        "export",
        "html",
        "--no-include-code",
        str(notebook_path),
        "-o",
        str(output_file),
        "--",
        str(data_file),
        str(links_file),
        "--fdr",
        str(fdr_threshold),
        "--genes",
        str(genes_mapped_threshold),
        "--max-terms",
        str(max_terms),
    ]

    logger.info(f"Running command: {shlex.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Marimo export completed: {output_file.name}")
        if result.stdout:
            logger.debug(f"Stdout:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"Stderr:\n{result.stderr}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Marimo export failed with return code {e.returncode}")
        if e.stdout:
            logger.error(f"Stdout:\n{e.stdout}")
        if e.stderr:
            logger.error(f"Stderr:\n{e.stderr}")
        return False


def prepare_data_input(data_dir: Path, output_dir: Path) -> Path:
    """
    Prepare data input by extracting zip files or copying directories.

    Args:
        data_dir: Input path (zip file or directory)
        output_dir: Output directory

    Returns:
        The actual data directory to use
    """
    if data_dir == output_dir:
        logger.info("Input and output directories are the same, no preparation needed")
        return data_dir

    if data_dir.is_file() and data_dir.suffix.lower() == ".zip":
        logger.info(f"Extracting zip file {data_dir} to {output_dir}")
        with zipfile.ZipFile(data_dir, "r") as zip_ref:
            zip_ref.extractall(output_dir)
        return output_dir

    elif data_dir.is_dir():
        logger.info(f"Copying directory {data_dir} to {output_dir}")
        if output_dir.exists():
            shutil.rmtree(output_dir)
        shutil.copytree(data_dir, output_dir)
        return output_dir

    else:
        raise ValueError(
            f"Data input must be either a zip file or directory: {data_dir}"
        )


def get_contrast_count(links_file: Path) -> int:
    """Get the number of contrasts from the links file."""
    with open(links_file) as f:
        return len(f.readlines())


def create_minimal_index(output_dir: Path) -> None:
    """Create a minimal index.html file with links to reports."""
    index_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>String GSEA Report (Marimo)</title>
    <style>
        body { font-family: system-ui, sans-serif; max-width: 800px; margin: 2rem auto; padding: 1rem; }
        h1 { color: #333; }
        ul { list-style-type: none; padding: 0; }
        li { margin: 0.5rem 0; }
        a { color: #0066cc; text-decoration: none; }
        a:hover { text-decoration: underline; }
    </style>
</head>
<body>
    <h1>STRING-GSEA Reports</h1>
    <ul>
        <li><a href="reports/index.html">Index - Overview and Links</a></li>
        <li><a href="reports/networks.html">Networks - Per-Contrast Visualizations</a></li>
        <li><a href="reports/multiple.html">Multiple Contrasts - Cross-Contrast Comparisons</a></li>
    </ul>
</body>
</html>"""
    index_path = output_dir / "index.html"
    index_path.write_text(index_content)
    logger.info(f"Created index file: {index_path}")


@app.default()
def render_marimo_reports(
    data_dir: Path,
    output_dir: Path | None = None,
    fdr_threshold: float = 0.05,
    genes_mapped_threshold: int = 10,
    max_terms: int = 100,
    zip: bool = False,
    reports_subfolder: bool = True,
) -> None:
    """
    Render Marimo reports as static HTML files.

    Args:
        data_dir: Path to data folder or zip file containing GSEA results.
        output_dir: Where the rendered HTML should go. Defaults to data_dir.
        fdr_threshold: FDR threshold for filtering results (default: 0.05).
        genes_mapped_threshold: Minimum genes mapped threshold (default: 10).
        max_terms: Maximum terms per category/contrast to render (default: 100).
        zip: Whether to create a zip archive of the output.
        reports_subfolder: Whether to put reports in a 'reports' subfolder.
    """
    data_dir = data_dir.absolute()
    if not data_dir.exists():
        raise FileNotFoundError(f"Data input does not exist: {data_dir}")

    if output_dir is None:
        output_dir = data_dir.with_suffix("") if data_dir.is_file() else data_dir
    else:
        output_dir = output_dir.absolute()

    if reports_subfolder:
        output_dir_dest = output_dir / "reports"
        output_dir_dest.mkdir(parents=True, exist_ok=True)
    else:
        output_dir_dest = output_dir
        output_dir_dest.mkdir(parents=True, exist_ok=True)

    # Prepare data input
    actual_data_dir = prepare_data_input(data_dir, output_dir_dest)

    # Find data files
    xlsx_files = sorted(actual_data_dir.glob("**/*_string_gsea_results_long.xlsx"))
    links_files = sorted(actual_data_dir.glob("**/*links.txt"))

    if not xlsx_files:
        raise FileNotFoundError(
            f"No XLSX files found matching '*_string_gsea_results_long.xlsx' in {actual_data_dir}"
        )
    if not links_files:
        raise FileNotFoundError(
            f"No links.txt files found in {actual_data_dir}"
        )

    if len(xlsx_files) != 1:
        logger.warning(f"Found {len(xlsx_files)} xlsx files, using first one")
    if len(links_files) != 1:
        logger.warning(f"Found {len(links_files)} links files, using first one")

    xlsx_file = xlsx_files[0]
    links_file = links_files[0]

    contrast_count = get_contrast_count(links_file)
    logger.info(f"Found {contrast_count} contrast(s)")

    # Get path to marimo reports
    with importlib.resources.path("string_gsea", "marimo_reports") as reports_path:
        # Export index report
        execute_marimo_export(
            notebook_path=reports_path / "report_index.py",
            output_file=output_dir_dest / "index.html",
            data_file=xlsx_file,
            links_file=links_file,
            fdr_threshold=fdr_threshold,
            genes_mapped_threshold=genes_mapped_threshold,
            max_terms=max_terms,
        )

        # Export networks report
        execute_marimo_export(
            notebook_path=reports_path / "report_networks.py",
            output_file=output_dir_dest / "networks.html",
            data_file=xlsx_file,
            links_file=links_file,
            fdr_threshold=fdr_threshold,
            genes_mapped_threshold=genes_mapped_threshold,
            max_terms=max_terms,
        )

        # Export multiple contrasts report (only if > 1 contrast)
        if contrast_count > 1:
            execute_marimo_export(
                notebook_path=reports_path / "report_multiple.py",
                output_file=output_dir_dest / "multiple.html",
                data_file=xlsx_file,
                links_file=links_file,
                fdr_threshold=fdr_threshold,
                genes_mapped_threshold=genes_mapped_threshold,
                max_terms=max_terms,
            )
        else:
            logger.info("Skipping multiple contrasts report (only 1 contrast)")

    if reports_subfolder:
        create_minimal_index(output_dir)

    if zip:
        zip_path = StringGSEAResults.zip_folder(output_dir)
        logger.info(f"Zipped reports to {zip_path}")

    logger.info(f"Reports rendered to {output_dir_dest}")


if __name__ == "__main__":
    app()
