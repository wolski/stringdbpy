import subprocess
from pathlib import Path
from cyclopts import App
from loguru import logger
import importlib.resources
import shlex
import shutil
import zipfile

from string_gsea.string_gsea_results import StringGSEAResults


app = App(
    help="Render Quarto reports in the package's docs/ directory, passing input_dir as a data_file param."
)


def execute_quarto_command(
    docs_path,
    output_dir,
    xlsx_file,
    links_file,
    FDR_threshold: float = 0.05,
    genes_mapped_threshold: int = 10,
    profile: str = "multiple",
):
    """
    Build and execute the quarto command with the given parameters.

    Args:
        docs_path: Path to the docs directory containing _quarto.yml
        output_dir: Directory where rendered HTML should go
        xlsx_file: Path to the xlsx results file
        links_file: Path to the links file
    """
    cmd = [
        "quarto",
        "render",
        str(docs_path),
        "--profile",
        profile,
        "--to",
        "html",
        "--output-dir",
        str(output_dir),
        "-P",
        f"data_file:{xlsx_file}",
        "-P",
        f"links_file:{links_file}",
        "-P",
        f"FDR_threshold:{FDR_threshold}",
        "-P",
        f"genes_mapped_threshold:{genes_mapped_threshold}",
    ]

    logger.info(f"Running command: {shlex.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("Quarto render completed successfully")
        if result.stdout:
            logger.debug(f"Stdout:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"Stderr:\n{result.stderr}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Quarto render failed with return code {e.returncode}")
        if e.stdout:
            logger.error(f"Stdout:\n{e.stdout}")
        if e.stderr:
            logger.error(f"Stderr:\n{e.stderr}")
        raise


def prepare_data_input(data_dir: Path, output_dir: Path) -> Path:
    """
    Prepare data input by extracting zip files or copying directories to output location.

    Args:
        data_dir: Input path that can be either a zip file or directory
        output_dir: Output directory where data should be placed

    Returns:
        Path: The actual data directory to use for finding xlsx and links files

    Raises:
        FileNotFoundError: If the input data_dir doesn't exist
        ValueError: If data_dir is neither a file nor directory
    """

    # If input and output are the same, no preparation needed
    if data_dir == output_dir:
        logger.info("Input and output directories are the same, no preparation needed")
        return data_dir

    if data_dir.is_file() and data_dir.suffix.lower() == ".zip":
        # Extract zip file to output directory
        logger.info(f"Extracting zip file {data_dir} to {output_dir}")
        with zipfile.ZipFile(data_dir, "r") as zip_ref:
            zip_ref.extractall(output_dir)
        return output_dir

    elif data_dir.is_dir():
        # Copy directory contents to output directory
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
    with open(links_file, "r") as f:
        return len(f.readlines())


# Copy required QMD files for single profile
def copy_qmd_files(docs_path: Path, output_dir_dest: Path, profile: str = "multiple"):
    single_qmd_files = [
        "_quarto.yml",
        "_quarto-single.yml",
        "index_single.qmd",
        "VisualizeNetworks.qmd",
    ]
    multiple_qmd_files = [
        "_quarto.yml",
        "_quarto-multiple.yml",
        "index_multiple.qmd",
        "VisualizeNetworks.qmd",
        "VisualizeMultipleContrastsGSEA.qmd",
    ]
    if profile == "single":
        qmd_files = single_qmd_files
    elif profile == "multiple":
        qmd_files = multiple_qmd_files
    else:
        raise ValueError(f"Invalid profile: {profile}")
    for filename in qmd_files:
        file_path = docs_path / filename
        if file_path.exists():
            shutil.copy2(file_path, output_dir_dest)


def create_minimal_index(output_dir: Path):
    """Create a minimal index.html file with link to reports."""
    index_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>String GSEA Report</title>
</head>
<body>
    <h1><a href="reports/index.html">String GSEA Report</a></h1>
</body>
</html>"""
    index_path = output_dir / "index.html"
    index_path.write_text(index_content)


@app.default()
def render_quarto_docs(
    data_dir: Path,
    output_dir: Path | None = None,
    FDR_threshold: float = 0.05,
    genes_mapped_threshold: int = 10,
    zip: bool = False,
    reports_subfolder: bool = True,
):
    """
    Render Quarto reports using data from the specified directory or zip file.

    Args:
        data_dir: path to your data folder or zip file (will be globbed for the .xlsx and links.txt).
        output_dir: where the rendered HTML should go. Defaults to data_dir/rendered_reports.
    """
    # Set up output directory
    data_dir = data_dir.absolute()
    if not data_dir.exists():
        raise FileNotFoundError(f"Data input does not exist: {data_dir}")

    if output_dir is None:
        output_dir = data_dir.with_suffix("") if data_dir.is_file() else data_dir
    else:
        output_dir = output_dir.absolute()

    if reports_subfolder:
        output_dir_dest = output_dir / "reports"  # / "rendered_reports"
        output_dir_dest.mkdir(parents=True, exist_ok=True)
    else:
        output_dir_dest = output_dir
        output_dir_dest.mkdir(parents=True, exist_ok=True)

    # Prepare data input (extract zip or copy directory)
    actual_data_dir = prepare_data_input(data_dir, output_dir_dest)

    # Find docs path and verify it contains required files
    with importlib.resources.path("string_gsea", "docs") as docs_path:
        xlsx_files = sorted(actual_data_dir.glob("**/*_string_gsea_results_long.xlsx"))
        links_files = sorted(actual_data_dir.glob("**/*links.txt"))

        # Log warnings if file counts are unexpected
        if len(xlsx_files) != 1:
            logger.warning(f"Found {len(xlsx_files)} xlsx files")
        if len(links_files) != 1:
            logger.warning(f"Found {len(links_files)} links files")

        contrast_count = get_contrast_count(links_files[0])
        profile = "multiple" if contrast_count > 1 else "single"

        copy_qmd_files(docs_path, output_dir_dest, profile)

        # Execute quarto command with the appropriate files
        execute_quarto_command(
            docs_path=output_dir_dest,
            output_dir=output_dir_dest,
            xlsx_file=xlsx_files[0],
            links_file=links_files[0],
            FDR_threshold=FDR_threshold,
            genes_mapped_threshold=genes_mapped_threshold,
            profile=profile,
        )

    if reports_subfolder:
        create_minimal_index(output_dir)

    if zip:
        # Zip the output directory
        zip_path = StringGSEAResults.zip_folder(output_dir)
        logger.info(f"Zipped reports to {zip_path}")


if __name__ == "__main__":
    app()
