import subprocess
from pathlib import Path
from cyclopts import App
from loguru import logger
import importlib.resources
import shlex
import shutil

from string_gsea.string_gsea_results import StringGSEAResults


app = App(
    help="Render Quarto reports in the package's docs/ directory, passing input_dir as a data_file param."
)

def execute_quarto_command(docs_path,
 output_dir,
 xlsx_file,
 links_file,
 FDR_threshold: float = 0.05,
 genes_mapped_threshold: int = 10):
    """
    Build and execute the quarto command with the given parameters.
    
    Args:
        docs_path: Path to the docs directory containing _quarto.yml
        output_dir: Directory where rendered HTML should go
        xlsx_file: Path to the xlsx results file
        links_file: Path to the links file
    """
    cmd = [
        "quarto", "render",
        str(docs_path),
        "--to", "html",
        "--output-dir", str(output_dir),
        "-P", f"data_file:{xlsx_file}",
        "-P", f"links_file:{links_file}",
        "-P", f"FDR_threshold:{FDR_threshold}",
        "-P", f"genes_mapped_threshold:{genes_mapped_threshold}"
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

@app.default()
def render_quarto_docs(data_dir: Path,
 output_dir: Path | None = None,
 FDR_threshold: float = 0.05,
 genes_mapped_threshold: int = 10,
 zip: bool = False
 ):
    """
    Render Quarto reports using data from the specified directory.
    
    Args:
        data_dir: path to your data folder (will be globbed for the .xlsx and links.txt).
        output_dir: where the rendered HTML should go. Defaults to data_dir/rendered_reports.
    """
    # Check if data_dir exists and make absolute path
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")
    data_dir = data_dir.absolute()

    # Set up output directory
    if output_dir is None:
        output_dir = data_dir 
    else:
        output_dir = output_dir.absolute()
    
    output_dir_dest = output_dir / "rendered_reports"
    output_dir_dest.mkdir(parents=True, exist_ok=True)

    # Find docs path and verify it contains required files
    with importlib.resources.path("string_gsea", "docs") as docs_path:
        if not (docs_path / "_quarto.yml").exists():
            raise FileNotFoundError(f"Missing _quarto.yml in: {docs_path}")
        # Copy .qmd and .yaml files to output directory
        for file in docs_path.glob("*.qmd"):
            shutil.copy2(file, output_dir_dest)
        for file in docs_path.glob("*.yml"):
            shutil.copy2(file, output_dir_dest)
        # Find data files
        xlsx_files = sorted(data_dir.glob("**/*_string_gsea_results_long.xlsx"))
        links_files = sorted(data_dir.glob("**/*links.txt"))

        # Log warnings if file counts are unexpected
        if len(xlsx_files) != 1:
            logger.warning(f"Found {len(xlsx_files)} xlsx files")
        if len(links_files) != 1:
            logger.warning(f"Found {len(links_files)} links files")

        # Execute quarto command with the appropriate files
        execute_quarto_command(
            docs_path=output_dir_dest,
            output_dir=output_dir_dest,
            xlsx_file=xlsx_files[0],
            links_file=links_files[0],
            FDR_threshold=FDR_threshold,
            genes_mapped_threshold=genes_mapped_threshold
        )

        if zip:
            # Zip the output directory
            zip_path = StringGSEAResults.zip_folder(output_dir)
            logger.info(f"Zipped reports to {zip_path}")

if __name__ == "__main__":
    app()