from cyclopts import App
from pathlib import Path
from loguru import logger

from string_gsea.gsea_utilities import find_zip_files
from string_gsea.run_string_gsea_bfabric import extract_workunit_id_from_file, run_string_gsea_bfabric

app = App()


@app.default()
def string_gsea_bfabric(
        in_dir: str,
        fdr: float = 0.25,
        out_dir: str = ".",
    ):
    """
    Run STRING GSEA analysis on the provided zip file.

    Args:
        in_dir (str): Path to the input zip file.
        fdr (float, optional): False discovery rate threshold. Defaults to 0.25.
        out_dir (str, optional): Base directory for output files. Defaults to ".".
    """

    in_dir = Path(in_dir)
    workunit_id = extract_workunit_id_from_file(in_dir / "params.yml")
    if workunit_id is None:
        raise ValueError("Workunit ID not found in params.yml")

    logger.info(f"Workunit ID: {workunit_id}")
    zip_path = find_zip_files(in_dir)[0]

    if not Path(zip_path).exists():
        raise FileNotFoundError(f"Zip file not found: {zip_path}")
    run_string_gsea_bfabric(zip_path, workunit_id, fdr=fdr, base_dir=Path(out_dir))

if __name__ == '__main__':
    app()
