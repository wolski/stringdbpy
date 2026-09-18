"""Polymorphic discovery of report template directories."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


@dataclass(frozen=True, slots=True)
class TemplatePaths:
    """Resolved Quarto template and vignette directories."""

    templates: Path
    vignettes: Path


class TemplateLocator(Protocol):
    """One source of installed or workspace report templates."""

    def locate(self) -> TemplatePaths | None:
        """Return usable paths, or ``None`` when this source is unavailable."""
        ...


class RPackageTemplateLocator:
    """Locate templates from an installed ``protsea`` R package."""

    @staticmethod
    def _system_file(directory: str) -> Path | None:
        if shutil.which("Rscript") is None:
            return None
        code = f"cat(system.file('{directory}', package='protsea'))"
        result = subprocess.run(
            ["Rscript", "-e", code],
            capture_output=True,
            text=True,
            check=False,
        )
        resolved = result.stdout.strip()
        return Path(resolved) if result.returncode == 0 and resolved else None

    def locate(self) -> TemplatePaths | None:
        templates = self._system_file("templates")
        if templates is None:
            return None
        return TemplatePaths(templates=templates, vignettes=templates)


@dataclass(frozen=True, slots=True)
class WorkspaceTemplateLocator:
    """Locate templates in the companion package checkout."""

    repository_root: Path

    def locate(self) -> TemplatePaths | None:
        package = self.repository_root.parent / "protsea"
        paths = TemplatePaths(
            templates=package / "inst" / "templates",
            vignettes=package / "inst" / "templates",
        )
        return paths if paths.templates.exists() and paths.vignettes.exists() else None


@dataclass(frozen=True, slots=True)
class OrderedTemplateLocator:
    """Select the first injected template source that is available."""

    locators: tuple[TemplateLocator, ...]

    def locate(self) -> TemplatePaths:
        for locator in self.locators:
            paths = locator.locate()
            if paths is not None:
                return paths
        raise FileNotFoundError("Could not locate protsea report templates")
