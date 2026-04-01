import zipfile
from pathlib import Path

from loguru import logger

from string_gsea.models.gsea_models import RankList, RankListCollection


def get_rank_files(zip_path) -> RankListCollection:
    """Read .rnk files from a ZIP archive and return a RankListCollection."""
    rank_lists = []
    with zipfile.ZipFile(zip_path, "r") as z:
        rnk_files = [f for f in z.namelist() if f.endswith(".rnk")]
        for file in rnk_files:
            with z.open(file) as f:
                content = f.read().decode()
            entries: dict[str, float] = {}
            for line in content.splitlines():
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                entries[parts[0]] = float(parts[1])
            filename_only = Path(file).stem
            rank_lists.append(RankList(contrast=filename_only, entries=entries))
    return RankListCollection(analysis="from_rnk", rank_lists=rank_lists)


def write_rank_files(ranks: RankListCollection, outdir: Path) -> Path:
    """Write .rnk files to disk, organized by analysis subdirectory."""
    sub = outdir / ranks.analysis
    sub.mkdir(parents=True, exist_ok=True)
    for rank_list in ranks:
        filepath = sub / f"{rank_list.contrast}.rnk"
        filepath.write_text(rank_list.to_rnk_string())
        logger.info(f"Wrote rank file: {filepath}")
    return outdir
