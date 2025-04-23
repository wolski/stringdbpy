import yaml
from dataclasses import dataclass, field
from pathlib import Path

@dataclass
class GSEASession:
    """
    Data class representing a STRING-GSEA session, with support for YAML serialization.
    Keys for res_job_id and res_data are tuple(str, str).
    """
    workunit_id: str
    species: int
    config_dict: dict
    base_path: Path
    res_job_id: dict[tuple[str, str], str] = field(default_factory=dict)
    res_data: dict[tuple[str, str], dict] = field(default_factory=dict)
    end_point_status = (
        "https://version-12-0.string-db.org/api/"
        "json/valuesranks_enrichment_status"
    )

    def to_yaml(self, filepath: Path = None) -> str:
        """
        Serialize this session to a YAML string or file, converting tuple keys to strings.
        """
        dump_dict = {
            'workunit_id': self.workunit_id,
            'species': self.species,
            'config_dict': self.config_dict,
            'base_path': str(self.base_path),
            'res_job_id': {f"{k[0]}~{k[1]}": v for k, v in self.res_job_id.items()},
            'res_data': {f"{k[0]}~{k[1]}": v for k, v in self.res_data.items()}
        }
        yaml_str = yaml.safe_dump(dump_dict, sort_keys=False)
        if filepath:
            with open(filepath, 'w') as f:
                f.write(yaml_str)
        return yaml_str

    @classmethod
    def from_yaml(cls, yaml_input) -> 'GSEASession':
        """
        Deserialize a session from a YAML string or file, converting string keys back to tuples.
        """
        if isinstance(yaml_input, (str, Path)) and Path(yaml_input).is_file():
            raw = Path(yaml_input).read_text()
        else:
            raw = yaml_input

        data = yaml.safe_load(raw)
        # Restore base_path
        data['base_path'] = Path(data['base_path'])
        # Convert res_job_id keys
        raw_jobs = data.get('res_job_id', {}) or {}
        data['res_job_id'] = {
            tuple(k.split('~', 1)): v for k, v in raw_jobs.items()
        }
        # Convert res_data keys
        raw_res = data.get('res_data', {}) or {}
        data['res_data'] = {
            tuple(k.split('~', 1)): v for k, v in raw_res.items()
        }
        return cls(
            workunit_id=data['workunit_id'],
            species=data['species'],
            config_dict=data['config_dict'],
            base_path=data['base_path'],
            res_job_id=data['res_job_id'],
            res_data=data['res_data']
        )
