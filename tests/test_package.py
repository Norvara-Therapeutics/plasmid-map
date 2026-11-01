"""Release invariants, checked here so a build never discovers them.

A version that disagrees with itself, an export that does not import, a module
missing from the distribution: each of these passes a source-tree run and
fails for whoever installs the package.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

import pytest

import plasmid_map

ROOT = Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "src" / "plasmid_map"


def configured() -> dict:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def _is_stdlib(name: str) -> bool:
    return name in sys.stdlib_module_names


def workflow(name: str) -> str:
    """One release workflow's text, or skip.

    The workflows are repository infrastructure and are deliberately not in
    the distribution, so these checks have nothing to read when the suite is
    run from an unpacked sdist. Skipping keeps that run green without
    weakening the check where it can be made.
    """
    path = ROOT / ".github" / "workflows" / name
    if not path.exists():
        pytest.skip("the release workflows are not shipped in the distribution")
    return path.read_text(encoding="utf-8")


class TestVersion:
    def test_the_two_version_literals_agree(self) -> None:
        """The tag, the metadata and the runtime all have to say one thing."""
        assert configured()["project"]["version"] == plasmid_map.__version__

    def test_the_version_is_a_prerelease_of_one_point_zero(self) -> None:
        assert plasmid_map.__version__ == "1.0.0rc2"


class TestExports:
    def test_every_exported_name_exists(self) -> None:
        missing = [
            name for name in plasmid_map.__all__ if not hasattr(plasmid_map, name)
        ]
        assert missing == []

    def test_the_export_list_is_sorted_and_unique(self) -> None:
        assert plasmid_map.__all__ == sorted(set(plasmid_map.__all__))

    def test_no_submodule_is_shadowed_by_an_export(self) -> None:
        """Binding a function over a module name breaks importing the module."""
        modules = {path.stem for path in PACKAGE.glob("*.py")} - {
            "__init__",
            "__main__",
        }
        shadowed = [name for name in plasmid_map.__all__ if name in modules]
        assert shadowed == []

    def test_the_two_draw_functions_are_reached_through_their_modules(self) -> None:
        """Both are called draw; the package cannot export either as that."""
        import plasmid_map.circular
        import plasmid_map.linear

        assert plasmid_map.circular.draw is not plasmid_map.linear.draw
        assert "draw" not in plasmid_map.__all__

    def test_the_tick_module_is_named_apart_from_the_function_it_exports(
        self,
    ) -> None:
        """``ruler`` is exported, so the module holding it is ``ticks``:
        binding the function over a module name breaks importing the module."""
        import plasmid_map.ticks

        assert plasmid_map.ticks.ruler is plasmid_map.ruler

    def test_every_module_can_be_imported_on_its_own(self) -> None:
        import importlib

        for path in sorted(PACKAGE.glob("*.py")):
            if path.stem == "__main__":
                continue
            importlib.import_module(f"plasmid_map.{path.stem}".rstrip("."))


class TestDistribution:
    def test_there_are_no_third_party_dependencies(self) -> None:
        assert configured()["project"]["dependencies"] == []

    def test_the_development_extra_can_run_what_the_readme_documents(self) -> None:
        """check_wheel.py builds, so the extra that installs it needs build."""
        extra = configured()["project"]["optional-dependencies"]["dev"]
        names = {item.split(">")[0].split("=")[0].split(";")[0].strip() for item in extra}
        assert {"pytest", "pyflakes", "build"} <= names

    def test_the_dev_extra_can_read_this_file_on_every_supported_python(self) -> None:
        """tomllib arrived in 3.11, and these tests read pyproject.toml.

        Dropping the marker fails at collection on 3.10, before any test
        runs, so nothing else here would catch it.
        """
        extra = configured()["project"]["optional-dependencies"]["dev"]
        assert any(
            item.startswith("tomli") and "python_version < '3.11'" in item
            for item in extra
        )

    def test_the_distribution_is_vendor_prefixed(self) -> None:
        assert configured()["project"]["name"] == "norvara-plasmid-map"

    def test_the_command_points_at_the_cli(self) -> None:
        scripts = configured()["project"]["scripts"]
        assert scripts == {"plasmid-map": "plasmid_map.cli:main"}

    def test_the_typing_marker_ships(self) -> None:
        assert (PACKAGE / "py.typed").exists()
        data = configured()["tool"]["setuptools"]["package-data"]
        assert data["plasmid_map"] == ["py.typed"]

    def test_the_manifest_carries_what_the_sdist_needs(self) -> None:
        """An unpacked sdist runs the suite, which reads the examples."""
        manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
        for line in (
            "include LICENSE",
            "include README.md",
            "recursive-include tests *.py",
            "recursive-include examples",
            "recursive-include scripts *.py",
        ):
            assert line in manifest
        for suffix in (".gb", ".tsv"):
            assert suffix in manifest

    def test_the_license_names_this_project(self) -> None:
        """A licence copied from a sibling repo keeps the sibling's notice."""
        licence = (ROOT / "LICENSE").read_text(encoding="utf-8")
        assert "Copyright 2026 plasmid-map contributors" in licence
        assert configured()["project"]["license"] == "PolyForm-Noncommercial-1.0.0"

    def test_the_python_floor_matches_the_tested_versions(self) -> None:
        assert configured()["project"]["requires-python"] == ">=3.10"

    def test_the_ci_matrix_starts_where_requires_python_does(self) -> None:
        """Testing only above the floor is how a 3.11-only import ships."""
        floor = configured()["project"]["requires-python"].lstrip(">=")
        found = re.search(r"python-version: \[([^\]]+)\]", workflow("prerelease.yml"))
        assert found is not None
        assert f'"{floor}"' in found.group(1)


class TestSources:
    @pytest.mark.parametrize(
        "path", sorted(PACKAGE.glob("*.py")), ids=lambda path: path.name
    )
    def test_every_module_is_documented(self, path: Path) -> None:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        assert ast.get_docstring(tree), f"{path.name} has no module docstring"

    @pytest.mark.parametrize(
        "path", sorted(PACKAGE.glob("*.py")), ids=lambda path: path.name
    )
    def test_every_public_function_is_documented(self, path: Path) -> None:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        undocumented = [
            node.name
            for node in tree.body
            if isinstance(node, (ast.FunctionDef, ast.ClassDef))
            and not node.name.startswith("_")
            and not ast.get_docstring(node)
        ]
        assert undocumented == []

    def test_nothing_imports_a_third_party_package(self) -> None:
        allowed = {"plasmid_map"}
        for path in PACKAGE.glob("*.py"):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for name in node.names:
                        root = name.name.split(".")[0]
                        assert _is_stdlib(root) or root in allowed, name.name
                elif isinstance(node, ast.ImportFrom) and node.module:
                    root = node.module.split(".")[0]
                    assert _is_stdlib(root) or root in allowed, node.module

    def test_the_coordinate_convention_is_stated_where_it_is_defined(self) -> None:
        """Every result depends on it, so the module that fixes it says so."""
        text = (PACKAGE / "feature.py").read_text(encoding="utf-8")
        assert "0-based and half-open" in text

    def test_the_table_convention_is_stated_where_the_tables_are_written(
        self,
    ) -> None:
        text = (PACKAGE / "rows.py").read_text(encoding="utf-8")
        assert "count from 1 and include both ends" in text


class TestScripts:
    def test_the_wheel_check_exists(self) -> None:
        assert (ROOT / "scripts" / "check_wheel.py").exists()

    def test_the_wheel_check_names_this_command(self) -> None:
        text = (ROOT / "scripts" / "check_wheel.py").read_text(encoding="utf-8")
        assert '"plasmid-map"' in text
        assert "bed-lite" not in text

    def test_the_example_builder_names_its_source(self) -> None:
        text = (ROOT / "examples" / "fetch_examples.py").read_text(encoding="utf-8")
        assert "eutils.ncbi.nlm.nih.gov" in text
        assert "--check" in text

    def test_the_example_checksums_are_recorded(self) -> None:
        """So that an example edited by hand is noticed."""
        sys.path.insert(0, str(ROOT / "examples"))
        import fetch_examples

        assert fetch_examples.RECORDS
        for name, value in fetch_examples.CHECKSUMS.items():
            assert len(value) == 64, f"no checksum recorded for {name}"

    def test_the_committed_examples_match_their_checksums(self) -> None:
        sys.path.insert(0, str(ROOT / "examples"))
        import fetch_examples

        assert fetch_examples.main(["--check"]) == 0

    def test_the_composed_example_is_declared_as_composed(self) -> None:
        """It is not a record and is not fetched, and says so."""
        sys.path.insert(0, str(ROOT / "examples"))
        import fetch_examples

        assert "cassette.tsv" in fetch_examples.COMPOSED
        assert "cassette.tsv" not in fetch_examples.RECORDS


class TestWorkflows:
    def test_the_release_workflows_name_this_package(self) -> None:
        for name in ("prerelease.yml", "publish-pypi.yml"):
            text = workflow(name)
            assert "src/plasmid_map/__init__.py" in text
            assert "plasmid_map" in text
            assert "bed_lite" not in text
            assert "align_pair" not in text

    def test_the_publish_default_tag_is_this_version(self) -> None:
        """A tag default copied from a sibling checks out a tag that is not here."""
        assert f"default: v{plasmid_map.__version__}" in workflow("publish-pypi.yml")

    def test_the_publisher_targets_this_distribution(self) -> None:
        assert "https://pypi.org/p/norvara-plasmid-map" in workflow(
            "publish-pypi.yml"
        )

    def test_the_source_check_covers_every_directory_of_python(self) -> None:
        for name in ("prerelease.yml", "publish-pypi.yml"):
            text = workflow(name)
            assert "src/plasmid_map tests scripts examples" in text
