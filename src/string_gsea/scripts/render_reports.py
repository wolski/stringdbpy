import subprocess
from pathlib import Path
from cyclopts import App, Positional

# Resolve package base directory
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = PACKAGE_ROOT / "docs"

def render_quarto_docs(data_dir: Path, output_dir: Path):
    if not DOCS_DIR.exists():
        raise FileNotFoundError(f"Docs directory does not exist: {DOCS_DIR}")
    if not (DOCS_DIR / "_quarto.yml").exists():
        raise FileNotFoundError(f"Missing _quarto.yml in: {DOCS_DIR}")

    # in the data_dir find all files with the suffix string_gsea_results_long.xlsx
    xlsx_files = list(data_dir.glob("*.xlsx"))



    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "quarto", "render",
        str(DOCS_DIR),
        "--to", "html",
        "--output-dir", str(output_dir),
        "--execute-params", f"data_file='{data_file}'"
    ]

    subprocess.run(cmd, check=True)


app = App(
    description="Render Quarto reports in the package's docs/ directory, passing input_dir as a data_file param."
)

@app.command
def run(
    input_dir: Path = Positional(help="Path to input data (used as param 'data_file' in QMDs)."),
    output_dir: Path = Positional(help="Path to output folder where HTML will be saved."),
):
    render_quarto_docs(input_dir, output_dir)


if __name__ == "__main__":
    app()
