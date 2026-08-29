"""Tests for template strategies, rendering, packaging, and command composition."""

import subprocess
import zipfile
from pathlib import Path

import pytest
import yaml

from string_gsea import workflow_cli
from string_gsea.workflow.packaging import package_results
from string_gsea.workflow.rendering import render_report
from string_gsea.workflow.templates import (
    OrderedTemplateLocator,
    RPackageTemplateLocator,
    TemplatePaths,
    WorkspaceTemplateLocator,
)

COMPATIBILITY = Path(__file__).parent / "data" / "compatibility"


def test_workspace_and_ordered_template_locators(tmp_path: Path) -> None:
    package = tmp_path / "stringGSEAplot"
    (package / "inst" / "templates").mkdir(parents=True)
    (package / "vignettes").mkdir()
    workspace = WorkspaceTemplateLocator(tmp_path)
    assert OrderedTemplateLocator((workspace,)).locate() == TemplatePaths(
        package / "inst" / "templates", package / "vignettes"
    )
    with pytest.raises(FileNotFoundError):
        OrderedTemplateLocator((WorkspaceTemplateLocator(tmp_path / "missing"),)).locate()


def test_r_package_locator_uses_return_codes(monkeypatch: pytest.MonkeyPatch) -> None:
    def executable(_name: str) -> str:
        return "/Rscript"

    monkeypatch.setattr("string_gsea.workflow.templates.shutil.which", executable)
    outputs = iter(("/templates\n", "/vignettes\n"))

    def run(*_args: object, **_kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess([], 0, stdout=next(outputs), stderr="")

    monkeypatch.setattr("string_gsea.workflow.templates.subprocess.run", run)
    assert RPackageTemplateLocator().locate() == TemplatePaths(
        Path("/templates"), Path("/vignettes")
    )


def test_rendering_copies_runs_and_cleans(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    dataset = tmp_path / "dataset"
    workunit = dataset / "WU_W1_GSEA"
    templates = tmp_path / "templates"
    vignettes = tmp_path / "vignettes"
    workunit.mkdir(parents=True)
    templates.mkdir()
    vignettes.mkdir()
    (vignettes / "GSEA_report.qmd").write_text("report")
    for name in ("index.qmd", "_fgcz-report.yml", "fgcz_header_quarto.html"):
        (templates / name).write_text(name)
    (workunit / "plots").mkdir()
    calls: list[list[str]] = []

    def run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("string_gsea.workflow.rendering.subprocess.run", run)
    done = workunit / "done.txt"
    render_report(dataset, "W1", TemplatePaths(templates, vignettes), done)
    assert len(calls) == 2
    assert done.exists()
    assert not (workunit / "GSEA_report.qmd").exists()
    assert not (workunit / "plots").exists()


def test_package_results_preserves_manifest_and_archive_shape(tmp_path: Path) -> None:
    workunit = tmp_path / "WU_W1_GSEA"
    workunit.mkdir()
    (workunit / "result.txt").write_text("result")
    manifest = tmp_path / "outputs.yml"
    archive = package_results(tmp_path, "W1", manifest)
    assert yaml.safe_load(manifest.read_text()) == {
        "outputs": [
            {
                "local_path": str(archive.resolve()),
                "store_entry_path": "WU_W1_GSEA.zip",
                "type": "bfabric_copy_resource",
            }
        ]
    }
    normalized_manifest = manifest.read_text().replace(str(tmp_path.resolve()), "<OUTPUT>")
    assert normalized_manifest == (COMPATIBILITY / "outputs.yml").read_text()
    with zipfile.ZipFile(archive) as zipped:
        assert zipped.namelist() == ["WU_W1_GSEA/", "WU_W1_GSEA/result.txt"]


def test_workflow_command_builds_dry_run_and_propagates_exit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    commands: list[list[str]] = []

    def run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        return subprocess.CompletedProcess(command, 7)

    paths = TemplatePaths(tmp_path / "templates", tmp_path / "vignettes")
    monkeypatch.setattr(workflow_cli, "_template_paths", lambda: paths)
    monkeypatch.setattr(workflow_cli.subprocess, "run", run)
    with pytest.raises(SystemExit, match="7"):
        workflow_cli.run_workflow("input.zip", "W1", "out", fdr=0.1, dry_run=True)
    assert commands[0][-1] == "-n"
    assert "fdr=0.1" in commands[0]
    assert f"templates_dir={paths.templates}" in commands[0]
    assert f"vignettes_dir={paths.vignettes}" in commands[0]
