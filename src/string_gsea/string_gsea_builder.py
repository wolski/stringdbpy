from datetime import datetime
from pathlib import Path

from string_gsea.gsea_config import GSEAConfig
from string_gsea.gsea_session import GSEASession
from string_gsea.models.gsea_models import RankListCollection
from string_gsea.string_api_client import StringAPIClient


class StringGSEABuilder:
    """
    Builder for submitting GSEA jobs to STRING-db and polling for completion.

    After ``submit().poll()``, access ``self.session`` for result data.
    """

    def __init__(
        self,
        rank_lists: RankListCollection,
        config: GSEAConfig,
        workunit_id: str = "12345",
        species: int = 9606,
        base_path: Path = Path("."),
    ):
        if not config:
            raise ValueError("config is None")

        self.session = GSEASession(
            current_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            workunit_id=workunit_id,
            species=species,
            config_dict=config,
            base_path=base_path,
            res_job_id={},
            res_data={},
        )
        self.rank_lists = rank_lists
        self.api_client = StringAPIClient(config, species)

    def submit(self) -> "StringGSEABuilder":
        analysis = self.rank_lists.analysis
        for rank_list in self.rank_lists:
            key = (analysis, rank_list.contrast)
            rank_str = rank_list.to_rnk_string()
            self.session.res_job_id[key] = self.api_client.submit_ranks(rank_str)
        return self

    def poll(self) -> "StringGSEABuilder":
        if not self.session.res_job_id:
            raise RuntimeError("No jobs to poll. Call submit() first.")
        for key, job_id in self.session.res_job_id.items():
            self.session.res_data[key] = self.api_client.poll_job(job_id)
        return self

    def save_session(self, filepath: Path) -> Path:
        filepath.parent.mkdir(parents=True, exist_ok=True)
        self.session.to_yaml(filepath)
        return filepath
