import subprocess
from pathlib import Path
from cyclopts import App
from loguru import logger

# Resolve package base directory
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = PACKAGE_ROOT / "docs"

app = App(
    help="Render Quarto reports in the package's docs/ directory, passing input_dir as a data_file param."
)
@app.default()
def render_quarto_docs(data_dir: Path, output_dir: Path):
    """
    input_dir: path to your data folder (will be globbed for the .xlsx and links.txt).
    output_dir: where the rendered HTML should go.
    """
    if not DOCS_DIR.exists():
        raise FileNotFoundError(f"Docs directory does not exist: {DOCS_DIR}")
    if not (DOCS_DIR / "_quarto.yml").exists():
        raise FileNotFoundError(f"Missing _quarto.yml in: {DOCS_DIR}")

    # in the data_dir find all files with the suffix string_gsea_results_long.xlsx
    xlsx_file = list([data_dir.glob("*_string_gsea_results_long.xlsx")]).sort()
    links_file = list([data_dir.glob("*links.txt")]).sort()

    if len(xlsx_file) != 1:
        logger.warning(f"Found {len(xlsx_file)} xlsx files and {len(links_file)} links files")
    if len(links_file) != 1:
        logger.warning(f"Found {len(links_file)} links files")

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "quarto", "render",
        str(DOCS_DIR),
        "--to", "html",
        "--output-dir", str(output_dir),
        "--execute-params", f"data_file='{xlsx_file[0]}' links_file='{links_file[0]}'"
    ]
    logger.info(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("Quarto render completed successfully")
        if result.stdout:
            logger.debug(f"Stdout:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"Stderr:\n{result.stderr}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Quarto render failed with return code {e.returncode}")
        if e.stdout:
            logger.error(f"Stdout:\n{e.stdout}")
        if e.stderr:
            logger.error(f"Stderr:\n{e.stderr}")
        raise





if __name__ == "__main__":
    app()
