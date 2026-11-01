"""The command line, run the way a person runs it."""

from __future__ import annotations

import json
import xml.etree.ElementTree as ElementTree
from pathlib import Path

import pytest

from plasmid_map.cli import build_parser, main


def run(argv, capsys) -> tuple[int, str, str]:
    """One command, with what it printed and what it complained about."""
    code = main(argv)
    captured = capsys.readouterr()
    return code, captured.out, captured.err


class TestCircular:
    def test_a_map_is_drawn_and_summarised(self, examples: Path, capsys) -> None:
        code, out, _err = run(["circular", str(examples / "pbr322.gb")], capsys)
        assert code == 0
        assert "drawn as a circular map" in out
        assert "4,361 bp" in out

    def test_the_svg_is_written_where_it_was_asked_for(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        target = tmp_path / "map.svg"
        code, _out, _err = run(
            ["circular", str(examples / "pbr322.gb"), "--svg", str(target)], capsys
        )
        assert code == 0
        body = target.read_text(encoding="utf-8")
        assert body.startswith("<?xml")
        ElementTree.fromstring(body.split("\n", 1)[1])

    def test_the_svg_can_go_to_standard_output(
        self, examples: Path, capsys
    ) -> None:
        code, out, _err = run(
            ["circular", str(examples / "pbr322.gb"), "-f", "svg"], capsys
        )
        assert code == 0
        assert out.startswith("<?xml")

    def test_a_report_is_a_whole_html_page(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        target = tmp_path / "map.html"
        run(
            [
                "circular",
                str(examples / "pbr322.gb"),
                "--report",
                str(target),
                "-f",
                "none",
            ],
            capsys,
        )
        assert target.read_text(encoding="utf-8").startswith("<!DOCTYPE html>")

    def test_the_features_are_written_as_a_table(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        target = tmp_path / "features.tsv"
        run(
            ["circular", str(examples / "pbr322.gb"), "-o", str(target), "-f", "none"],
            capsys,
        )
        lines = target.read_text(encoding="utf-8").splitlines()
        assert lines[0].split("\t")[0] == "kind"
        assert len(lines) == 49

    def test_json_carries_the_summary_and_the_notes(
        self, examples: Path, capsys
    ) -> None:
        code, out, _err = run(
            ["circular", str(examples / "pbr322.gb"), "-f", "json"], capsys
        )
        assert code == 0
        payload = json.loads(out)
        assert payload["construct"]["length"] == 4361
        assert payload["summary"]["features_drawn"] == 48
        assert payload["notes"]

    def test_nothing_is_printed_when_nothing_is_asked_for(
        self, examples: Path, capsys
    ) -> None:
        code, out, _err = run(
            ["circular", str(examples / "pbr322.gb"), "-f", "none"], capsys
        )
        assert code == 0 and out == ""

    def test_the_map_can_be_turned(self, examples: Path, capsys) -> None:
        first = main(["circular", str(examples / "pbr322.gb"), "-f", "svg"])
        turned = capsys.readouterr().out
        assert first == 0
        main(
            ["circular", str(examples / "pbr322.gb"), "-f", "svg", "--start", "2000"]
        )
        assert capsys.readouterr().out != turned


class TestLinear:
    def test_a_track_is_drawn(self, examples: Path, capsys) -> None:
        code, out, _err = run(["linear", str(examples / "aav2.gb")], capsys)
        assert code == 0
        assert "drawn as a linear map" in out

    def test_turning_a_linear_track_is_not_an_option_it_takes(
        self, examples: Path, capsys
    ) -> None:
        """There is nothing to turn on a line, so it is refused rather than
        offered and quietly ignored."""
        with pytest.raises(SystemExit) as raised:
            main(["linear", str(examples / "aav2.gb"), "--start", "100"])
        assert raised.value.code == 2
        assert "--start" in capsys.readouterr().err

    def test_a_circular_record_is_cut_and_the_cut_is_reported(
        self, examples: Path, capsys
    ) -> None:
        code, out, _err = run(
            ["linear", str(examples / "phix174.gb"), "--kind", "gene"], capsys
        )
        assert code == 0
        assert "two pieces" in out


class TestWithoutAFile:
    def test_a_table_and_a_length_are_enough(
        self, examples: Path, capsys
    ) -> None:
        code, out, _err = run(
            [
                "linear",
                "--length",
                "2744",
                "--features",
                str(examples / "cassette.tsv"),
            ],
            capsys,
        )
        assert code == 0
        assert "9 in 2 rows" in out

    def test_annotations_alone_are_enough(self, capsys) -> None:
        code, out, _err = run(
            [
                "circular",
                "--length",
                "5000",
                "--annotate",
                "EGFP:1021..1740:+:CDS",
                "--annotate",
                "ori:3000..3600",
                "-f",
                "json",
            ],
            capsys,
        )
        assert code == 0
        assert json.loads(out)["summary"]["features_drawn"] == 2

    def test_a_length_with_nothing_on_it_is_refused(self, capsys) -> None:
        code, _out, err = run(["circular", "--length", "5000"], capsys)
        assert code == 2
        assert "--annotate" in err

    def test_neither_a_file_nor_a_length_is_refused(self, capsys) -> None:
        code, _out, err = run(["circular"], capsys)
        assert code == 2
        assert "--length" in err

    def test_a_length_with_a_file_is_refused(
        self, examples: Path, capsys
    ) -> None:
        code, _out, err = run(
            ["circular", str(examples / "pbr322.gb"), "--length", "100"], capsys
        )
        assert code == 2
        assert "states its own length" in err

    def test_a_gc_track_without_a_sequence_is_refused(self, capsys) -> None:
        code, _out, err = run(
            ["circular", "--length", "5000", "--annotate", "x:1..100", "--gc"], capsys
        )
        assert code == 2
        assert "carries none" in err


class TestSelection:
    def test_a_kind_filter_narrows_the_map(self, examples: Path, capsys) -> None:
        code, out, _err = run(
            ["circular", str(examples / "pbr322.gb"), "--kind", "CDS", "-f", "json"],
            capsys,
        )
        assert code == 0
        assert json.loads(out)["summary"]["features_drawn"] == 3

    def test_hiding_a_kind_leaves_the_rest(
        self, examples: Path, capsys
    ) -> None:
        code, out, _err = run(
            ["circular", str(examples / "pbr322.gb"), "--hide", "misc_binding"], capsys
        )
        assert code == 0
        assert "misc_binding" not in out.split("note:")[0]

    def test_a_minimum_length_leaves_the_short_ones_out(
        self, examples: Path, capsys
    ) -> None:
        code, out, _err = run(
            [
                "circular",
                str(examples / "pbr322.gb"),
                "--minimum-length",
                "100",
            ],
            capsys,
        )
        assert code == 0
        assert "--minimum-length" in out

    def test_filtering_everything_away_is_refused(
        self, examples: Path, capsys
    ) -> None:
        code, _out, err = run(
            ["circular", str(examples / "pbr322.gb"), "--kind", "tRNA"], capsys
        )
        assert code == 2
        assert "nothing left to draw" in err


class TestFeaturesCommand:
    def test_it_lists_what_would_be_drawn(self, examples: Path, capsys) -> None:
        code, out, _err = run(["features", str(examples / "pbr322.gb")], capsys)
        assert code == 0
        assert "would be drawn" in out

    def test_it_can_write_the_whole_table(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        target = tmp_path / "f.tsv"
        run(
            ["features", str(examples / "pbr322.gb"), "-o", str(target)], capsys
        )
        assert len(target.read_text(encoding="utf-8").splitlines()) == 49

    def test_json_is_a_structure(self, examples: Path, capsys) -> None:
        code, out, _err = run(
            ["features", str(examples / "pbr322.gb"), "-f", "json"], capsys
        )
        assert code == 0
        assert json.loads(out)["length"] == 4361

    def test_a_long_table_is_cut_short_for_a_terminal(
        self, examples: Path, capsys
    ) -> None:
        code, out, _err = run(["features", str(examples / "pbr322.gb")], capsys)
        assert code == 0
        assert "and 23 more" in out


class TestTypesCommand:
    def test_it_lists_the_kinds_and_their_colours(self, capsys) -> None:
        code, out, _err = run(["types"], capsys)
        assert code == 0
        assert "CDS" in out and "#3d6fb4" in out

    def test_the_mono_palette_can_be_listed(self, capsys) -> None:
        code, out, _err = run(["types", "--palette", "mono"], capsys)
        assert code == 0
        assert "mono palette" in out

    def test_an_override_shows_in_the_listing(self, capsys) -> None:
        code, out, _err = run(["types", "--colour", "CDS=#123456"], capsys)
        assert code == 0
        assert "#123456" in out

    def test_json_is_a_structure(self, capsys) -> None:
        code, out, _err = run(["types", "-f", "json"], capsys)
        assert code == 0
        assert json.loads(out)["palette"] == "colour"


class TestOutputPaths:
    def test_writing_over_the_input_is_refused(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        source = tmp_path / "x.gb"
        source.write_text(
            (examples / "pbr322.gb").read_text(encoding="utf-8"), encoding="utf-8"
        )
        code, _out, err = run(
            ["circular", str(source), "--svg", str(source)], capsys
        )
        assert code == 2
        assert "being read" in err

    def test_two_outputs_to_one_path_are_refused(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        target = tmp_path / "both"
        code, _out, err = run(
            [
                "circular",
                str(examples / "pbr322.gb"),
                "--svg",
                str(target),
                "-o",
                str(target),
            ],
            capsys,
        )
        assert code == 2
        assert "different paths" in err

    def test_the_features_command_will_not_overwrite_its_input(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        source = tmp_path / "x.gb"
        source.write_text(
            (examples / "pbr322.gb").read_text(encoding="utf-8"), encoding="utf-8"
        )
        code, _out, err = run(["features", str(source), "-o", str(source)], capsys)
        assert code == 2
        assert "being read" in err

    def test_writing_over_the_feature_table_is_refused(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        table = tmp_path / "f.tsv"
        table.write_text("name\tstart\tend\nx\t10\t20\n", encoding="utf-8")
        code, _out, err = run(
            [
                "linear",
                "--length",
                "100",
                "--features",
                str(table),
                "-o",
                str(table),
            ],
            capsys,
        )
        assert code == 2
        assert "feature table" in err


class TestErrors:
    def test_a_file_that_is_not_there_is_reported_not_raised(
        self, capsys
    ) -> None:
        code, _out, err = run(["circular", "/no/such/file.gb"], capsys)
        assert code == 2
        assert err.startswith("plasmid-map: ")

    def test_a_colour_that_is_not_one_is_refused_by_the_parser(self) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(["types", "--colour", "CDS=blurple"])

    def test_a_colour_override_needs_a_kind_and_a_colour(self) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(["types", "--colour", "nonsense"])

    def test_a_width_of_zero_is_refused_by_the_parser(self) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(["circular", "x.gb", "--width", "0"])

    def test_a_map_too_narrow_for_its_labels_is_refused_with_advice(
        self, examples: Path, capsys
    ) -> None:
        code, _out, err = run(
            ["circular", str(examples / "pbr322.gb"), "--width", "150"], capsys
        )
        assert code == 2
        assert "--label-width" in err


class TestParser:
    def test_the_two_drawing_commands_take_the_same_options(self) -> None:
        parser = build_parser()
        commands = parser._subparsers._group_actions[0].choices
        circular = {
            flag
            for action in commands["circular"]._actions
            for flag in action.option_strings
        }
        linear = {
            flag
            for action in commands["linear"]._actions
            for flag in action.option_strings
        }
        assert circular - linear == {"--start"}
        assert linear - circular == set()

    def test_every_command_is_offered(self) -> None:
        parser = build_parser()
        commands = parser._subparsers._group_actions[0].choices
        assert list(commands) == ["circular", "linear", "features", "types"]

    def test_the_version_is_the_package_version(self, capsys) -> None:
        from plasmid_map import __version__

        with pytest.raises(SystemExit):
            main(["--version"])
        assert __version__ in capsys.readouterr().out
