import os
import shutil
import tempfile
import requests
import time
from loguru import logger
from pathlib import Path
from datetime import datetime
from string_gsea.get_species import OxFieldsZip
from string_gsea.gsea_session import GSEASession
from string_gsea.gsea_utilities import get_rank_files
from string_gsea.string_gsea_results import StringGSEAResults
from string_gsea.gsea_config import GSEAConfig
from dataclasses import asdict

class StringGSEABuilder:
    """
    Builder for submitting GSEA jobs to STRING-db, polling for completion,
    constructing a results object, and serializing/deserializing the build state via YAML.
    """

    def __init__(
        self,
        rank_dataframes: dict,
        config: GSEAConfig,
        workunit_id: str = "12345",
        species: int = 9606,
        base_path: Path = Path('.')
    ):
        if not config:
            raise ValueError("config is None")

        # Initialize session data class
        self.session = GSEASession(
            current_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            workunit_id=workunit_id,
            species=species,
            config_dict=config,
            base_path=base_path,
            res_job_id={},
            res_data={}
        )
        # Store rank inputs separately
        self.rank_dataframes = rank_dataframes

    def _submit_single(self, rank_data: str) -> str:
        base = "https://version-12-0.string-db.org/api"
        url = f"{base}/json/valuesranks_enrichment_submit"
        params = {
            "species": self.session.species,
            "caller_identity": asdict(self.session.config_dict)['caller_identity'],
            "identifiers": rank_data,
            "api_key": asdict(self.session.config_dict)['api_key'],
            "ge_fdr": asdict(self.session.config_dict)['fdr'],
            "ge_enrichment_rank_direction" : asdict(self.session.config_dict)['ge_enrichment_rank_direction'],
        }
        resp = requests.post(url, data=params)
        resp.raise_for_status()
        data = resp.json()[0]
        if data.get('status') == 'error':
            msg = data.get('message', 'no message')
            logger.error(f"STRING-db error: {msg}")
            raise RuntimeError(f"STRING-db error: {msg}")
        job_id = data['job_id']
        logger.info(f"Job submitted: {job_id}")
        return job_id

    def submit(self) -> 'StringGSEABuilder':
        for key, df in self.rank_dataframes.items():
            rank_str = df.write_csv(separator='\t', include_header=False)
            self.session.res_job_id[key] = self._submit_single(rank_str)
        return self

    def _poll_single(self, job_id: str, sleep_t: int = 10, max_time: int = 3600) -> dict:
        status_url = (
            self.session.end_point_status
            + f"?api_key={asdict(self.session.config_dict)['api_key']}&job_id={job_id}"
        )
        elapsed = 0
        while elapsed < max_time:
            resp = requests.get(status_url)
            resp.raise_for_status()
            record = resp.json()[0]
            st = record.get('status')
            logger.info(f"Polling {job_id}: status={st} after {elapsed}s")
            if st == 'success':
                return record
            if st in ('nothing found', 'unknown organism'):
                raise RuntimeError(f"Polling error {st} for job {job_id}")
            time.sleep(sleep_t)
            elapsed += sleep_t
        raise TimeoutError(f"Job {job_id} did not finish in {max_time}s")

    def poll(self) -> 'StringGSEABuilder':
        if not self.session.res_job_id:
            raise RuntimeError("No jobs to poll. Call submit() first.")
        for key, job_id in self.session.res_job_id.items():
            self.session.res_data[key] = self._poll_single(job_id)
        return self

    def write_rank_files(self) -> Path:
        outdir = self.get_res_path()
        for (outer, inner), df in self.rank_dataframes.items():
            sub = outdir / outer
            sub.mkdir(parents=True, exist_ok=True)
            filepath = sub / f"{inner}.rnk"
            with open(filepath, 'wb') as f:
                f.write(df.write_csv(separator='\t', include_header=False).encode())
            logger.info(f"Wrote rank file: {filepath}")
        return outdir

    def get_res_path(self) -> Path:
        p = self.session.base_path / f"WU_{self.session.workunit_id}_GSEA"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def build_results(self) -> StringGSEAResults:
        """
        Wrap the current session into a StringGSEAResults for writing outputs.
        """
        # GSEASession.res_data must already be populated (via poll())
        return StringGSEAResults(self.session)

    def get_result(self) -> StringGSEAResults:
        return self.submit().poll().build_results()

    def save_session(self, filepath: Path = None) -> Path:
        out_dir = self.get_res_path()
        path = filepath or (out_dir / 'gsea_session.yml')
        self.session.to_yaml(path)
        return path


if __name__ == '__main__':
    config = GSEAConfig(
        api_key="b36F8oaRJwFZ",
        fdr=0.25,
        caller_identity="www.fgcz.ch",
        ge_enrichment_rank_direction=-1
    )
    print("Current working directory:", os.getcwd())
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    zip_path = project_root / "tests" / "data" / "2848501.zip"
    workunit_id = "abcd"
    species = OxFieldsZip.get_species_from_oxes(zip_path)

    logger.info(f"species is : {species}")
    dataframes = get_rank_files(zip_path)
    tempdir = Path(tempfile.mkdtemp())
    logger.info(f"Created temporary directory: {tempdir}")

    builder = StringGSEABuilder(
        rank_dataframes=dataframes,
        config=config,
        workunit_id=workunit_id,
        species=species,
        base_path=tempdir
    )
    builder.write_rank_files()
    builder.submit().poll()
    session_path = builder.save_session()
    results = builder.build_results()

    logger.info(f"Jobs: {builder.session.res_job_id}")
    results.write_links()
    results.write_gsea_tsv()
    results.write_gsea_graphs()
    builder.save_session()
    #results.save_session()
    # copy session_path file into tests/data/dummy_d
    shutil.copy(session_path, project_root / "tests" / "data" / "dummy_d" / "session.yml")
    results2 = StringGSEAResults(GSEASession.from_yaml(session_path))
    results2.write_links()


    print("Workflow complete.")
