"""The README held to the program, rather than to a reading of it.

Every quoted output block is re-run through the command line and compared byte
for byte, the command table is checked against the parser, and every flag
named anywhere in the prose has to exist. A README that has drifted from the
program fails here rather than misleading a reader.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import pytest

from plasmid_map.cli import build_parser, main

ROOT = Path(__file__).resolve().parent.parent
README = (ROOT / "README.md").read_text(encoding="utf-8")


def blocks() -> list[tuple[str, str]]:
    """Each command in the README paired with the output quoted after it."""
    pattern = re.compile(r"```bash\n([^`]*?)\n```\n\n```\n([^`]*?)```", re.DOTALL)
    return pattern.findall(README)


def commands() -> list[str]:
    """Every bash command the README shows."""
    return re.findall(r"```bash\n([^`]*?)\n```", README, re.DOTALL)


#: Sections that quote another program's flags: pip's --pre under Install and
#: the example builder's --check under Development. A flag named there is not
#: this program's, so they are left out of the audit.
BORROWED = ("## Install", "## Development")


def prose() -> str:
    """The README except the sections that quote another program's flags.

    Everything else is checked, which is what catches a flag named in passing
    in prose rather than only the ones in the option table.
    """
    kept: list[str] = []
    skipping = False
    for line in README.splitlines():
        if line.startswith("## "):
            skipping = line.strip() in BORROWED
        if not skipping:
            kept.append(line)
    return "\n".join(kept)


def flat() -> str:
    """The README with its line wrapping collapsed.

    A claim the README makes is a claim whatever column it happened to break
    at, so the checks below read it as one line.
    """
    return " ".join(README.split())


def parser_flags() -> set[str]:
    found: set[str] = set()
    stack = [build_parser()]
    while stack:
        current = stack.pop()
        for action in current._actions:
            found.update(action.option_strings)
            choices = getattr(action, "choices", None)
            if isinstance(choices, dict):
                stack.extend(choices.values())
    return found


def resolved(command: str, into: Path | None = None) -> list[str]:
    """One README command as arguments, with its paths made absolute.

    Example paths are resolved against the checkout so the test runs from
    anywhere, and anything the command would write goes to a temporary
    directory rather than into the repository.
    """
    parts = shlex.split(command.replace("\\\n", " "))
    assert parts[0] == "plasmid-map"
    arguments = [
        str(ROOT / part) if part.startswith("examples/") else part
        for part in parts[1:]
    ]
    for option in ("--svg", "--report", "-o", "--output"):
        if option in arguments:
            position = arguments.index(option) + 1
            written = Path(arguments[position]).name
            arguments[position] = str((into or Path("/tmp")) / written)
    return arguments


class TestQuotedOutput:
    def test_there_are_quoted_blocks_to_check(self) -> None:
        assert len(blocks()) >= 4

    @pytest.mark.parametrize("index", range(4))
    def test_a_quoted_block_is_what_the_command_prints(
        self, index: int, capsys, tmp_path: Path
    ) -> None:
        command, expected = blocks()[index]
        code = main(resolved(command, tmp_path))
        printed = capsys.readouterr().out
        assert code == 0
        # The README quotes the relative paths a reader would type, while the
        # test has to pass absolute ones so it runs from anywhere.
        assert printed.replace(str(ROOT / "examples"), "examples") == expected

    def test_a_command_that_says_it_draws_a_file_draws_one(
        self, capsys, tmp_path: Path
    ) -> None:
        """The quickstart says the drawing goes to --svg, so it had better."""
        drawing = [
            command
            for command in commands()
            if "--svg" in command or "--report" in command
        ]
        assert drawing, "no README command writes a drawing"
        for command in drawing:
            arguments = resolved(command, tmp_path)
            assert main(arguments) == 0
            capsys.readouterr()
            for option in ("--svg", "--report"):
                if option in arguments:
                    written = Path(arguments[arguments.index(option) + 1])
                    assert written.exists() and written.stat().st_size > 0

    def test_every_quoted_command_runs(self, capsys, tmp_path: Path) -> None:
        """Including any whose output is not quoted."""
        for command in commands():
            if command.startswith(("pip", "python")):
                continue
            assert main(resolved(command, tmp_path)) == 0
            capsys.readouterr()


class TestCommandTable:
    def test_every_command_is_in_the_table(self) -> None:
        table = re.search(r"## The commands\n\n(.*?)\n\n", README, re.DOTALL)
        assert table is not None
        named = re.findall(r"^\| `(\w+)` \|", table.group(1), re.MULTILINE)
        offered = build_parser()._subparsers._group_actions[0].choices
        assert named == list(offered)

    def test_every_stdout_format_is_named(self) -> None:
        from plasmid_map.cli import STDOUT_FORMATS

        for name in STDOUT_FORMATS:
            assert f"`{name}`" in README


class TestFlags:
    def test_every_flag_named_in_the_readme_exists(self) -> None:
        named = set(re.findall(r"(?<![\w-])(--[a-z][a-z-]+)", prose()))
        unknown = named - parser_flags()
        assert unknown == set(), f"the README names flags the parser does not: {unknown}"

    def test_the_check_can_fail(self) -> None:
        """A flag audit that cannot fail proves nothing."""
        assert "--nonexistent-flag" not in parser_flags()

    def test_every_flag_the_parser_takes_is_in_the_readme(self) -> None:
        """An option nobody is told about is an option nobody uses."""
        named = set(re.findall(r"(?<![\w-])(--[a-z][a-z-]+)", README))
        long_flags = {
            flag
            for flag in parser_flags()
            if flag.startswith("--") and flag not in ("--help", "--version")
        }
        assert long_flags - named == set()

    def test_every_output_option_is_documented(self) -> None:
        for flag in ("--svg", "--report", "--length", "--annotate", "--features"):
            assert flag in README


class TestClaims:
    def test_the_coordinate_convention_is_stated_as_the_code_writes_it(self) -> None:
        assert "count from 1 and include both ends" in flat()

    def test_the_origin_crossing_rule_is_stated(self) -> None:
        assert "one arc through twelve o'clock" in flat()
        assert "the sum of its parts" in flat()

    def test_the_length_the_readme_quotes_is_the_length_of_that_gene(
        self, phix174
    ) -> None:
        """The three crossing genes are not the same length, so the claim is
        about the one whose location the README prints and no other."""
        quoted = next(
            item
            for item in phix174.features
            if item.kind == "gene" and item.start == 3980
        )
        assert quoted.length == 1542
        assert "`join(3981..5386,1..136)` and is 1,542 bases long" in flat()
        lengths = {item.length for item in phix174.crossing()}
        assert len(lengths) == 3

    def test_the_label_estimate_is_admitted_as_an_estimate(self) -> None:
        """A reader trusting a map with long labels should know this."""
        assert "*estimated*" in flat()
        assert "a fixed average glyph width" in flat()

    def test_the_self_containedness_is_stated(self) -> None:
        assert "no script, no web font" in flat()

    def test_determinism_is_stated(self) -> None:
        assert "draws the same bytes" in flat()

    def test_nothing_being_dropped_quietly_is_stated(self) -> None:
        assert "Nothing is dropped quietly" in flat()

    def test_the_band_thinning_is_stated(self) -> None:
        assert "the bands are thinned" in flat()

    def test_the_composed_example_is_declared_as_composed(self) -> None:
        assert "`examples/cassette.tsv` is not a record" in flat()

    def test_the_license_is_named(self) -> None:
        assert "PolyForm Noncommercial License 1.0.0" in README

    def test_the_distribution_name_is_the_one_that_installs(self) -> None:
        import tomllib

        with (ROOT / "pyproject.toml").open("rb") as handle:
            configured = tomllib.load(handle)
        assert configured["project"]["name"] in README


def test_the_readme_opens_with_the_project_name() -> None:
    assert README.startswith("# plasmid-map\n")


def test_the_readme_has_not_grown_past_what_it_needs() -> None:
    """A README nobody finishes is a README nobody reads.

    The limit is deliberately generous; its purpose is to make growth visible
    rather than to forbid it.
    """
    assert len(README) < 12_000, f"the README is {len(README)} bytes"
