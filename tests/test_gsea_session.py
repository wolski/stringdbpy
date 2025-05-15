import pytest
import yaml
from pathlib import Path

from string_gsea.gsea_session import GSEASession


@pytest.fixture
def sample_session(tmp_path):
    # Create a sample GSEASession with tuple keys
    cfg = {
        'current_date': '2025-05-15 10:00:00',
        'api_key': 'test_key',
        'fdr': 0.05,
        'caller_identity': 'pytest',
        'ge_enrichment_rank_direction': 1
    }
    session = GSEASession(
        current_date='2025-05-15 10:00:00',
        workunit_id='WU_test',
        species=9606,
        config_dict=cfg,
        base_path=tmp_path,
        res_job_id={('outer', 'inner'): 'job123'},
        res_data={('outer', 'inner'): {'status': 'success', 'page_url': 'http://example.com'}}
    )
    return session


def test_to_yaml_string(sample_session):
    yaml_str = sample_session.to_yaml()
    data = yaml.safe_load(yaml_str)

    # Verify core fields
    assert data['workunit_id'] == sample_session.workunit_id
    assert data['species'] == sample_session.species
    assert data['config_dict'] == sample_session.config_dict
    # Serialized keys use 'outer~inner'
    assert data['res_job_id'] == {'outer~inner': 'job123'}
    assert data['res_data'] == {'outer~inner': {'status': 'success', 'page_url': 'http://example.com'}}


def test_to_yaml_file(tmp_path, sample_session):
    yaml_file = tmp_path / 'session.yml'
    yaml_str = sample_session.to_yaml(filepath=yaml_file)

    # File should exist and content should match returned string
    assert yaml_file.exists()
    file_text = yaml_file.read_text()
    assert file_text == yaml_str


def test_from_yaml_string(sample_session):
    yaml_str = sample_session.to_yaml()
    loaded = GSEASession.from_yaml(yaml_str)

    # Dataclass equality should hold
    assert loaded == sample_session


def test_from_yaml_file(tmp_path, sample_session):
    yaml_file = tmp_path / 'session.yml'
    sample_session.to_yaml(filepath=yaml_file)
    loaded = GSEASession.from_yaml(yaml_input=yaml_file)
    assert loaded == sample_session


def test_endpoint_status():
    # Ensure the class attribute is set correctly
    assert hasattr(GSEASession, 'end_point_status')
    assert isinstance(GSEASession.end_point_status, str)
    assert 'valuesranks_enrichment_status' in GSEASession.end_point_status
