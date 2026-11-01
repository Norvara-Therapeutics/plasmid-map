"""One test for each defect found in review, named for what it was.

These are kept apart from the tests of each module so that the record of what
was once wrong, and is now checked, stays readable. Every one of them failed
before the fix beside it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from conftest import construct, feature, sequence

from plasmid_map.cli import main
from plasmid_map.composition import measure
from plasmid_map.drawing import Style, groups_of, select
from plasmid_map.feature import Construct, Feature, Span, location_text, within
from plasmid_map.genbank import parse, parse_locus
from plasmid_map.lanes import assign
from plasmid_map.linear import draw as draw_linear
from plasmid_map.palette import Palette
from plasmid_map.rows import FEATURE_COLUMNS, feature_columns, feature_rows, kind_rows
from plasmid_map.tables import TableError, looks_like_gff3, read, read_gff3
from plasmid_map.ticks import coordinates


class TestSilentlyMissingTrack:
    def test_asking_for_skew_on_a_linear_map_draws_one(self) -> None:
        """It validated the profile, kept it, and then drew only GC."""
        held = construct(
            topology="linear",
            features=(feature(kind="CDS", start=0, end=500),),
            sequence=sequence(),
        )
        drawing = draw_linear(
            held, select(held), Style(skew=True), measure(held.sequence, False)
        )
        assert drawing.svg.count('class="skew-up"') > 0


class TestFalseOriginCrossing:
    def test_a_widened_terminal_feature_does_not_appear_at_the_start(self) -> None:
        """Widening a one base feature at the last base pushed its span past
        the end, and a span past the end was read as crossing the origin."""
        held = construct(
            length=5000,
            topology="linear",
            features=(feature(kind="misc_feature", start=4999, end=5000),),
        )
        drawing = draw_linear(held, select(held), Style())
        assert drawing.svg.count('class="feature"') == 1


class TestCutDirection:
    def test_the_point_is_on_the_end_the_feature_ends_at(self) -> None:
        held = construct(
            features=(feature(name="wrap", start=950, end=1050, strand=1),)
        )
        drawing = draw_linear(held, select(held), Style(labels="none"))
        import re

        shapes = re.findall(r'<path d="(M[^"]+)"[^>]*class="feature"', drawing.svg)
        pointed = [shape for shape in shapes if shape.count("L") == 4]
        assert len(pointed) == 1
        assert float(pointed[0][1:].split()[0]) < drawing.width / 2

    def test_no_join_is_drawn_through_the_middle_of_the_sequence(self) -> None:
        held = construct(features=(feature(name="wrap", start=950, end=1050),))
        assert 'class="join"' not in draw_linear(
            held, select(held), Style()
        ).svg


class TestGff3:
    def test_rows_sharing_an_id_are_one_feature(self) -> None:
        text = (
            "##gff-version 3\n"
            "p1\tx\tCDS\t1\t100\t.\t+\t0\tID=cds1;Name=rep\n"
            "p1\tx\tCDS\t200\t300\t.\t+\t0\tID=cds1;Name=rep\n"
        )
        features = read_gff3(text, 1000, False)
        assert len(features) == 1
        assert len(features[0].parts) == 2
        assert features[0].length == 201

    def test_rows_without_an_id_stay_separate(self) -> None:
        text = (
            "##gff-version 3\n"
            "p1\tx\tCDS\t1\t100\t.\t+\t0\tName=a\n"
            "p1\tx\tCDS\t200\t300\t.\t+\t0\tName=b\n"
        )
        assert len(read_gff3(text, 1000, False)) == 2

    def test_two_sequences_are_refused_rather_than_added_together(self) -> None:
        text = (
            "##gff-version 3\n"
            "chr1\tx\tCDS\t1\t100\t.\t+\t0\tName=a\n"
            "chr2\tx\tCDS\t200\t300\t.\t+\t0\tName=b\n"
        )
        with pytest.raises(TableError) as raised:
            read_gff3(text, 1000, False)
        assert "--features-seqid" in str(raised.value)

    def test_naming_one_sequence_keeps_only_its_rows(self) -> None:
        text = (
            "##gff-version 3\n"
            "chr1\tx\tCDS\t1\t100\t.\t+\t0\tName=a\n"
            "chr2\tx\tCDS\t200\t300\t.\t+\t0\tName=b\n"
        )
        features = read_gff3(text, 1000, False, seqid="chr2")
        assert [item.name for item in features] == ["b"]

    def test_a_sequence_the_file_does_not_name_is_refused(self) -> None:
        text = "##gff-version 3\nchr1\tx\tCDS\t1\t100\t.\t+\t0\tName=a\n"
        with pytest.raises(TableError) as raised:
            read_gff3(text, 1000, False, seqid="chr9")
        assert "chr1" in str(raised.value)

    def test_eight_columns_are_refused(self) -> None:
        with pytest.raises(TableError) as raised:
            read_gff3("##gff-version 3\np1\tx\tCDS\t1\t10\t.\t+\t0\n", 100, False)
        assert "9 tab-separated" in str(raised.value)

    def test_ten_columns_are_refused(self) -> None:
        with pytest.raises(TableError):
            read_gff3(
                "##gff-version 3\np1\tx\tCDS\t1\t10\t.\t+\t0\tName=a\textra\n",
                100,
                False,
            )

    def test_a_nine_column_table_with_a_header_is_not_gff3(self) -> None:
        """Nine columns is not enough to be GFF3; the columns have to read
        like GFF3 columns."""
        header = "name\tkind\tstart\tend\tstrand\ta\tb\tc\td\n"
        assert not looks_like_gff3(header + "x\tCDS\t1\t10\t+\t1\t2\t3\t4\n")
        features = read(header + "x\tCDS\t1\t10\t+\t1\t2\t3\t4\n", 100, False)
        assert features[0].name == "x"


class TestTopology:
    def test_a_locus_line_with_no_topology_is_not_a_statement(self) -> None:
        assert not parse_locus("LOCUS       X    100 bp    DNA").stated
        assert parse_locus(
            "LOCUS       X    100 bp    DNA     circular"
        ).stated

    def test_a_record_that_does_not_say_is_reported_as_not_saying(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "quiet.gb"
        path.write_text(
            "LOCUS       QUIET    40 bp    DNA\n"
            "FEATURES             Location/Qualifiers\n"
            "     CDS             1..30\n"
            '                     /gene="x"\n'
            "//\n",
            encoding="utf-8",
        )
        from plasmid_map.reader import load

        loaded = load(str(path))
        assert not loaded.topology_stated

    def test_making_a_record_circular_joins_what_meets_at_the_origin(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "linear.gb"
        path.write_text(
            "LOCUS       LIN    100 bp    DNA     linear\n"
            "FEATURES             Location/Qualifiers\n"
            "     CDS             join(90..100,1..10)\n"
            '                     /gene="wrap"\n'
            "//\n",
            encoding="utf-8",
        )
        from plasmid_map.reader import load

        linear_first = load(str(path)).first
        assert len(linear_first.features[0].parts) == 2
        made_circular = load(str(path), topology="circular").first
        assert len(made_circular.features[0].parts) == 1
        assert made_circular.features[0].crosses_origin(100)


class TestColours:
    def test_a_legend_tells_an_overridden_kind_apart_from_its_group(self) -> None:
        palette = Palette("colour", {"CDS": "#ff0000"})
        layout = assign(
            [
                feature(kind="CDS", name="a", start=0, end=100),
                feature(kind="gene", name="b", start=200, end=300),
            ],
            1000,
            True,
            1,
        )
        legend = dict(groups_of(layout, palette))
        assert legend["CDS"] == "#ff0000"
        assert legend["coding"] != "#ff0000"

    def test_types_lists_a_kind_that_only_exists_as_an_override(self) -> None:
        rows = kind_rows(Palette("colour", {"custom": "#ff0000"}))
        assert ["custom", "other", "#ff0000", 0] in [row for row in rows]

    def test_an_overridden_kind_is_not_called_a_default(self) -> None:
        from plasmid_map.text import notes

        held = construct(features=(feature(kind="custom", start=0, end=100),))
        drawing = draw_linear(
            held,
            select(held),
            Style(palette=Palette("colour", {"custom": "#ff0000"})),
        )
        assert not any("default" in note for note in notes(drawing))


class TestTables:
    def test_a_site_between_bases_does_not_end_before_it_starts(self) -> None:
        """Which is what this table uses to mean crossing the origin."""
        layout = assign(
            [Feature(kind="misc_feature", name="s", parts=(Span(102, 102),))],
            4361,
            True,
            1,
        )
        row = feature_rows(layout.placements, 4361)[0]
        assert (row[2], row[3]) == (102, 103)
        assert row[FEATURE_COLUMNS.index("length")] == 0
        assert row[FEATURE_COLUMNS.index("location")] == "102^103"

    def test_which_end_is_partial_survives_into_the_table(self) -> None:
        held = Feature(
            kind="CDS", name="p", parts=(Span(0, 300),), partial=(True, False)
        )
        layout = assign([held], 1000, False, 1)
        row = feature_rows(layout.placements, 1000)[0]
        assert row[FEATURE_COLUMNS.index("partial_ends")] == "start"
        assert row[FEATURE_COLUMNS.index("location")] == "<1..300"

    def test_a_complement_writes_its_marker_where_the_record_had_it(self) -> None:
        from plasmid_map.genbank import parse_location

        where = parse_location("complement(<1..300)")
        held = Feature(
            kind="CDS",
            name="p",
            parts=where.parts,
            strand=where.strand,
            partial=where.partial,
        )
        assert location_text(held, 1000) == "complement(<1..300)"

    def test_a_qualifier_can_be_carried_into_a_column(self) -> None:
        held = Feature(
            kind="CDS",
            name="x",
            parts=(Span(0, 100),),
            qualifiers=(("product", "a protein"),),
        )
        layout = assign([held], 1000, False, 1)
        columns = feature_columns(["product"])
        row = feature_rows(layout.placements, 1000, ["product"])[0]
        assert columns[-1] == "product"
        assert row[-1] == "a protein"


class TestGenBank:
    def test_a_repeated_qualifier_names_the_feature_the_way_it_is_read(
        self,
    ) -> None:
        """The label took the last; the qualifier API returned the first."""
        text = (
            "LOCUS       X    100 bp    DNA     circular\n"
            "FEATURES             Location/Qualifiers\n"
            "     misc_feature    1..10\n"
            '                     /note="first"\n'
            '                     /note="second"\n'
            "//\n"
        )
        held = parse(text)[0].features[0]
        assert held.name == held.qualifier("note") == "first"

    def test_a_wrapped_translation_holds_no_spaces(self) -> None:
        text = (
            "LOCUS       X    100 bp    DNA     circular\n"
            "FEATURES             Location/Qualifiers\n"
            "     CDS             1..30\n"
            '                     /translation="MKWVTFISLLLLFSSAYSRGVFRRDTHKS\n'
            '                     EIAHRFKDLGE"\n'
            "//\n"
        )
        held = parse(text)[0].features[0]
        assert " " not in (held.qualifier("translation") or " ")

    def test_a_site_between_the_last_base_and_the_first_is_a_place(self) -> None:
        held = Feature(kind="misc_feature", name="o", parts=(Span(100, 100),))
        assert within(held, 100, circular=True)
        assert not within(held, 100, circular=False)
        assert location_text(held, 100) == "100^1"

    def test_a_real_record_translation_is_amino_acids_only(self, pbr322) -> None:
        for item in pbr322.features:
            translation = item.qualifier("translation")
            if translation:
                assert translation.isalpha(), item.name


class TestTicks:
    def test_a_one_base_ruler_marks_base_one_once(self) -> None:
        assert coordinates(1, 1, circular=False) == [1]

    def test_a_short_ruler_does_not_repeat_itself(self) -> None:
        marked = coordinates(5, 1, circular=False)
        assert marked == sorted(set(marked))


class TestCommandLine:
    def test_two_things_cannot_be_written_into_one_stream(
        self, examples: Path, capsys
    ) -> None:
        code = main(["circular", str(examples / "pbr322.gb"), "--svg", "-", "-f", "json"])
        assert code == 2
        assert "one stream" in capsys.readouterr().err

    def test_one_thing_to_standard_output_is_fine(
        self, examples: Path, capsys
    ) -> None:
        code = main(["circular", str(examples / "pbr322.gb"), "--svg", "-", "-f", "none"])
        assert code == 0
        assert capsys.readouterr().out.startswith("<?xml")

    def test_an_option_that_would_do_nothing_is_refused(
        self, examples: Path, capsys
    ) -> None:
        code = main(
            ["circular", str(examples / "pbr322.gb"), "--gc-window", "200", "-f", "none"]
        )
        assert code == 2
        assert "--gc-window" in capsys.readouterr().err

    def test_a_feature_format_without_a_feature_table_is_refused(
        self, examples: Path, capsys
    ) -> None:
        code = main(
            [
                "circular",
                str(examples / "pbr322.gb"),
                "--features-format",
                "gff3",
                "-f",
                "none",
            ]
        )
        assert code == 2
        assert "--features-format" in capsys.readouterr().err

    def test_an_input_format_with_no_input_is_refused(self, capsys) -> None:
        code = main(
            [
                "circular",
                "--length",
                "1000",
                "--annotate",
                "x:1..100",
                "--input-format",
                "fasta",
                "-f",
                "none",
            ]
        )
        assert code == 2
        assert "--input-format" in capsys.readouterr().err

    def test_a_record_of_a_multi_record_file_can_be_chosen(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        both = tmp_path / "two.gb"
        both.write_text(
            (examples / "pbr322.gb").read_text(encoding="utf-8")
            + (examples / "phix174.gb").read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        assert main(["circular", str(both), "--record", "2", "-f", "json"]) == 0
        payload = json.loads(capsys.readouterr().out)
        assert payload["construct"]["length"] == 5386

    def test_a_record_can_be_chosen_by_accession(
        self, examples: Path, tmp_path: Path, capsys
    ) -> None:
        both = tmp_path / "two.gb"
        both.write_text(
            (examples / "pbr322.gb").read_text(encoding="utf-8")
            + (examples / "phix174.gb").read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        assert main(["circular", str(both), "--record", "NC_001422", "-f", "json"]) == 0
        assert json.loads(capsys.readouterr().out)["construct"]["length"] == 5386

    def test_a_record_that_is_not_there_is_refused_with_what_is(
        self, examples: Path, capsys
    ) -> None:
        code = main(["circular", str(examples / "pbr322.gb"), "--record", "9"])
        assert code == 2
        assert "holding 1" in capsys.readouterr().err

    def test_added_features_are_reported_in_every_format(
        self, examples: Path, capsys
    ) -> None:
        """Not only in the one a person reads."""
        assert (
            main(
                [
                    "linear",
                    "--length",
                    "2744",
                    "--features",
                    str(examples / "cassette.tsv"),
                    "-f",
                    "json",
                ]
            )
            == 0
        )
        payload = json.loads(capsys.readouterr().out)
        assert payload["summary"]["features_added"] == 9
        assert any("added to the record" in note for note in payload["notes"])

    def test_a_qualifier_column_reaches_the_table(
        self, examples: Path, capsys
    ) -> None:
        assert (
            main(
                [
                    "features",
                    str(examples / "pbr322.gb"),
                    "--kind",
                    "CDS",
                    "--qualifiers",
                    "product",
                    "-f",
                    "tsv",
                ]
            )
            == 0
        )
        printed = capsys.readouterr().out.splitlines()
        assert printed[0].split("\t")[-1] == "product"
        assert any(line.split("\t")[-1] for line in printed[1:])


class TestComposition:
    def test_a_circular_record_drawn_as_a_line_is_measured_as_a_line(
        self, examples: Path, capsys
    ) -> None:
        """Measuring round the origin and drawing on a line put three windows
        on the same rightmost point."""
        assert (
            main(
                [
                    "linear",
                    str(examples / "phix174.gb"),
                    "--kind",
                    "gene",
                    "--gc",
                    "-f",
                    "json",
                ]
            )
            == 0
        )
        payload = json.loads(capsys.readouterr().out)
        assert payload["summary"]["gc_measured_round_the_origin"] is False
        record = parse(
            (examples / "phix174.gb").read_text(encoding="utf-8")
        )[0]
        as_drawn = measure(
            record.sequence,
            False,
            payload["summary"]["gc_window"],
            payload["summary"]["gc_step"],
        )
        assert payload["summary"]["gc_windows"] == len(as_drawn.windows)
        assert as_drawn.windows[-1].end == record.length

    def test_the_summary_says_how_the_windows_were_measured(self) -> None:
        from plasmid_map.text import measures

        held = construct(features=(feature(start=0, end=100),), sequence=sequence())
        circular = draw_linear(
            held, select(held), Style(gc=True), measure(held.sequence, True)
        )
        linear = draw_linear(
            held, select(held), Style(gc=True), measure(held.sequence, False)
        )
        assert "round the origin" in dict(measures(circular))["windows"]
        assert "stopping at the end" in dict(measures(linear))["windows"]

    def test_a_skew_only_profile_is_not_called_gc_content(self) -> None:
        from plasmid_map.text import measures

        held = construct(features=(feature(start=0, end=100),), sequence=sequence())
        drawing = draw_linear(
            held, select(held), Style(skew=True), measure(held.sequence, False)
        )
        labels = dict(measures(drawing))
        assert "GC content" not in labels
        assert "GC skew" in labels


class TestGrammar:
    def test_one_feature_is_not_are(self) -> None:
        from plasmid_map.text import notes

        held = construct(length=5000, features=(feature(start=100, end=101),))
        drawing = draw_linear(held, select(held), Style())
        told = " ".join(notes(drawing))
        assert "1 feature are" not in told
        assert "1 feature is" in told

    def test_one_added_feature_is_not_were(self, capsys) -> None:
        code = main(
            ["circular", "--length", "5000", "--annotate", "cut:1200"]
        )
        assert code == 0
        printed = capsys.readouterr().out
        assert "1 feature were" not in printed
        assert "1 feature was added" in printed


class TestLabels:
    def test_a_label_wider_than_the_track_is_cut_to_fit(self) -> None:
        from plasmid_map.labels import Along, stack_linear

        rows = stack_linear(
            [Along(text="a very long feature name indeed", x=50.0)],
            top=0.0,
            font_size=11.0,
            left=0.0,
            right=100.0,
        )
        assert rows.labels[0].box.right <= 100.0 + 1e-9
        assert rows.shortened == 1

    def test_the_cutting_is_reported(self) -> None:
        from plasmid_map.text import notes

        held = Construct(
            name="t",
            length=1000,
            topology="linear",
            features=(
                Feature(
                    kind="CDS",
                    name="an extremely long feature name that will not fit",
                    parts=(Span(0, 1000),),
                ),
            ),
        )
        drawing = draw_linear(held, select(held), Style(width=120))
        assert any("cut further" in note for note in notes(drawing))
