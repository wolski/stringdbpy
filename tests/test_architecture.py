# pyright: reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnknownArgumentType=false
"""Executable directed-folder architecture invariants."""

from pathlib import Path

import grimp


def test_feature_packages_have_no_direct_sibling_dependencies() -> None:
    graph = grimp.build_graph("string_gsea")
    siblings = {
        "string_gsea.gsea",
        "string_gsea.ora",
        "string_gsea.taxonomy",
        "string_gsea.workflow",
    }
    for sibling in siblings:
        modules = {
            module
            for module in graph.modules
            if module == sibling or module.startswith(f"{sibling}.")
        }
        imported = set().union(
            *(graph.find_modules_directly_imported_by(module) for module in modules)
        )
        targets = {
            candidate
            for candidate in siblings - {sibling}
            if any(
                dependency == candidate or dependency.startswith(f"{candidate}.")
                for dependency in imported
            )
        }
        assert len(targets) <= 1, (sibling, targets)


def test_package_import_graph_has_no_cycles() -> None:
    graph = grimp.build_graph("string_gsea")
    assert graph.nominate_cycle_breakers("string_gsea") == set()


def test_package_initializers_are_empty() -> None:
    source = Path(__file__).parents[1] / "src" / "string_gsea"
    for initializer in source.rglob("__init__.py"):
        assert not initializer.read_text().strip(), initializer


def test_packaged_workflow_imports_only_current_modules() -> None:
    snakefile = (
        Path(__file__).parents[1] / "src" / "string_gsea" / "workflow" / "Snakefile"
    ).read_text()
    assert "string_gsea.workflow.helpers" not in snakefile
    assert "from string_gsea" not in snakefile
