"""The summary and the notes beside it."""

from __future__ import annotations

from conftest import construct, feature, sequence, spread

from plasmid_map.circular import draw
from plasmid_map.composition import measure
from plasmid_map.drawing import Style, select
from plasmid_map.linear import draw as draw_linear
from plasmid_map.palette import Palette
from plasmid_map.text import (
    bases,
    cell,
    definitions,
    measures,
    notes,
    number,
    palette_table,
    plural,
    share,
    summary,
    table,
    truncate,
)


def _drawing(**style):
    held = construct(features=spread(), name="pMAP1")
    return draw(held, select(held), Style(**style))


class TestNumbers:
    def test_a_share_that_rounds_to_zero_is_not_shown_as_zero(self) -> None:
        """Nothing found and almost nothing found say opposite things."""
        assert share(0.0001) == "<0.1%"

    def test_a_share_that_rounds_to_everything_is_not_shown_as_everything(
        self,
    ) -> None:
        assert share(0.99999) == ">99.9%"

    def test_an_ordinary_share_is_a_percentage(self) -> None:
        assert share(0.5374) == "53.7%"

    def test_a_share_that_does_not_exist_is_a_word(self) -> None:
        assert share(None) == "n/a"

    def test_a_whole_number_keeps_its_separators(self) -> None:
        assert number(4361.0) == "4,361"

    def test_a_count_and_its_noun_agree(self) -> None:
        assert plural(1, "label") == "1 label"
        assert plural(2, "label") == "2 labels"

    def test_a_length_is_abbreviated_once_it_needs_to_be(self) -> None:
        assert bases(999) == "999 bp"
        assert bases(4361) == "4.4 kb"
        assert bases(4_000_000) == "4.00 Mb"


class TestCells:
    def test_a_value_that_does_not_exist_is_blank_not_the_word_none(self) -> None:
        assert cell(None) == ""

    def test_a_boolean_is_a_word_a_person_reads(self) -> None:
        assert cell(True) == "yes"

    def test_a_float_is_rounded_for_a_terminal(self) -> None:
        assert cell(1 / 3) == "0.333"

    def test_a_long_value_is_cut(self) -> None:
        assert cell("x" * 80).endswith("\u2026")

    def test_whitespace_is_collapsed(self) -> None:
        assert truncate("a\nb  c") == "a b c"


class TestTable:
    def test_numeric_columns_are_right_aligned(self) -> None:
        written = table(["name", "count"], [["a", 1], ["bb", 100]], numeric=(1,))
        lines = written.splitlines()
        assert lines[2].endswith("1")
        assert lines[3].endswith("100")

    def test_a_table_with_no_rows_is_empty(self) -> None:
        assert table(["a"], []) == ""

    def test_definitions_line_their_labels_up(self) -> None:
        written = definitions([("a", "1"), ("bbb", "2")])
        assert written.splitlines()[0].startswith("a  ")


class TestMeasures:
    def test_the_length_of_the_record_is_reported(self) -> None:
        labels = dict(measures(_drawing()))
        assert "features in the record" in labels

    def test_the_lanes_are_named_rings_on_a_circle(self) -> None:
        assert "lane" in dict(measures(_drawing()))["drawn"]

    def test_the_lanes_are_named_rows_on_a_line(self) -> None:
        held = construct(topology="linear", features=spread()[:-1])
        drawing = draw_linear(held, select(held), Style())
        assert "row" in dict(measures(drawing))["drawn"]

    def test_the_gc_figures_appear_when_a_track_was_drawn(self) -> None:
        held = construct(features=spread(), sequence=sequence())
        drawing = draw(
            held, select(held), Style(gc=True, width=900), measure(held.sequence, True)
        )
        assert "GC content" in dict(measures(drawing))


class TestNotes:
    def test_the_whole_record_feature_being_left_out_is_said(self, pbr322) -> None:
        drawing = draw(pbr322, select(pbr322), Style(width=900))
        assert any("whole record" in note for note in notes(drawing))

    def test_features_left_out_by_a_filter_are_counted_by_kind(
        self, pbr322
    ) -> None:
        drawing = draw(
            pbr322, select(pbr322, kinds=["CDS"]), Style(width=900)
        )
        told = " ".join(notes(drawing))
        assert "kind filter" in told
        assert "misc_binding" in told

    def test_crossing_the_origin_is_said_and_the_features_named(
        self, phix174
    ) -> None:
        drawing = draw(
            phix174, select(phix174, kinds=["gene"]), Style(width=900)
        )
        told = " ".join(notes(drawing))
        assert "cross the origin" in told
        assert "one arc" in told

    def test_a_linear_track_says_the_crossing_features_were_cut(
        self, phix174
    ) -> None:
        drawing = draw_linear(
            phix174, select(phix174, kinds=["gene"]), Style(width=900)
        )
        assert any("two pieces" in note for note in notes(drawing))

    def test_widening_is_said_with_the_width_used(self) -> None:
        held = construct(length=5000, features=(feature(start=100, end=101),))
        drawing = draw(held, select(held), Style())
        assert any("drawn that wide" in note for note in notes(drawing))

    def test_moved_labels_are_counted(self) -> None:
        crowded = tuple(
            feature(name=f"feature {number}", start=number, end=number + 1)
            for number in range(20)
        )
        drawing = draw(
            construct(features=crowded), select(construct(features=crowded)), Style()
        )
        assert any("moved" in note for note in notes(drawing))

    def test_an_unstated_topology_is_said(self) -> None:
        drawing = _drawing()
        assert any("does not record" in note for note in notes(drawing, stated=False))

    def test_further_records_in_the_file_are_said(self) -> None:
        drawing = _drawing()
        assert any("records" in note for note in notes(drawing, loaded_records=3))

    def test_a_kind_with_no_colour_of_its_own_is_said(self) -> None:
        held = construct(features=(feature(kind="something_new", start=0, end=100),))
        drawing = draw(held, select(held), Style())
        assert any("default" in note for note in notes(drawing))

    def test_thinned_bands_are_said(self) -> None:
        many = tuple(
            feature(name=f"f{number}", start=0, end=900) for number in range(12)
        )
        drawing = draw(
            construct(features=many), select(construct(features=many)), Style(width=600)
        )
        assert any("thinned" in note for note in notes(drawing))

    def test_a_quiet_map_says_little(self) -> None:
        held = construct(
            topology="linear", features=(feature(kind="CDS", start=100, end=500),)
        )
        drawing = draw_linear(held, select(held), Style())
        assert notes(drawing) == []


class TestSummary:
    def test_the_heading_names_the_source_and_the_shape(self) -> None:
        written = summary(_drawing(), "x.gb", "genbank")
        assert written.startswith("x.gb (genbank, circular, 1,000 bp)")
        assert "circular map" in written.splitlines()[0]

    def test_every_note_is_marked_as_one(self) -> None:
        written = summary(_drawing(), "x.gb", "genbank")
        for note in notes(_drawing()):
            assert f"note: {note}" in written

    def test_the_kinds_drawn_are_tabled(self) -> None:
        assert "kind" in summary(_drawing(), "x.gb", "genbank")

    def test_it_ends_with_one_newline(self) -> None:
        written = summary(_drawing(), "x.gb", "genbank")
        assert written.endswith("\n") and not written.endswith("\n\n")


class TestPaletteTable:
    def test_every_kind_is_listed(self) -> None:
        from plasmid_map.palette import GROUPS

        written = palette_table(Palette())
        assert written.count("\n") >= len(GROUPS)
