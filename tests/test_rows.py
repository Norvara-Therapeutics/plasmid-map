"""The tables, and the coordinate convention they are written in."""

from __future__ import annotations

from conftest import construct, feature, sequence, spread

from plasmid_map.circular import draw
from plasmid_map.composition import measure
from plasmid_map.drawing import Style, select
from plasmid_map.feature import minimum_width
from plasmid_map.lanes import assign
from plasmid_map.rows import (
    FEATURE_COLUMNS,
    KIND_COLUMNS,
    WINDOW_COLUMNS,
    coordinate,
    counts,
    feature_rows,
    field,
    kind_rows,
    location,
    render,
    summary_rows,
    window_rows,
)
from plasmid_map.palette import Palette


def _rows(features, length=1000, circular=True):
    layout = assign(features, length, circular, minimum_width(length))
    return feature_rows(layout.placements, length)


class TestFields:
    def test_nothing_is_an_empty_field(self) -> None:
        assert field(None) == ""

    def test_a_float_keeps_its_precision(self) -> None:
        assert field(1 / 3) == repr(1 / 3)

    def test_a_whole_float_is_written_as_a_whole_number(self) -> None:
        assert field(5.0) == "5"

    def test_a_boolean_is_a_word(self) -> None:
        assert field(True) == "true" and field(False) == "false"

    def test_a_tab_inside_a_name_becomes_a_space(self) -> None:
        """So that one row stays one row."""
        assert "\t" not in field("a\tname")

    def test_a_newline_inside_a_name_becomes_a_space(self) -> None:
        assert "\n" not in field("a\nname")


class TestRender:
    def test_a_table_has_a_header_and_a_trailing_newline(self) -> None:
        written = render(("a", "b"), [[1, 2]])
        assert written == "a\tb\n1\t2\n"

    def test_a_table_with_no_rows_is_still_a_header(self) -> None:
        assert render(("a", "b"), []) == "a\tb\n"


class TestCoordinates:
    def test_the_first_base_is_one(self) -> None:
        assert coordinate(0, 1000) == 1

    def test_a_position_past_the_end_comes_round(self) -> None:
        assert coordinate(1050, 1000) == 51

    def test_a_feature_is_written_as_the_record_wrote_it(self) -> None:
        rows = _rows([feature(kind="CDS", name="bla", start=85, end=1276)], 4361)
        assert rows[0][2] == 86
        assert rows[0][3] == 1276

    def test_a_feature_reaching_the_last_base_says_so(self) -> None:
        rows = _rows([feature(start=0, end=1000)], 1000)
        assert (rows[0][2], rows[0][3]) == (1, 1000)

    def test_an_origin_crossing_feature_ends_before_it_starts(self) -> None:
        """As it does in a GenBank file, with a column that says so."""
        rows = _rows([feature(name="wrap", start=950, end=1050)])
        assert rows[0][2] == 951
        assert rows[0][3] == 50
        assert rows[0][FEATURE_COLUMNS.index("crosses_origin")] is True

    def test_its_length_is_the_bases_it_holds(self) -> None:
        rows = _rows([feature(name="wrap", start=950, end=1050)])
        assert rows[0][FEATURE_COLUMNS.index("length")] == 100


class TestLocation:
    def test_a_simple_range(self) -> None:
        assert location(feature(start=85, end=1276), 4361) == "86..1276"

    def test_a_single_base(self) -> None:
        assert location(feature(start=425, end=426), 4361) == "426"

    def test_a_site_between_two_bases(self) -> None:
        assert location(feature(start=102, end=102), 4361) == "102^103"

    def test_the_minus_strand(self) -> None:
        assert location(
            feature(start=3292, end=4153, strand=-1), 4361
        ) == "complement(3293..4153)"

    def test_a_spliced_feature(self) -> None:
        assert location(
            feature(parts=((320, 1906), (2227, 2252))), 4679
        ) == "join(321..1906,2228..2252)"

    def test_an_origin_crossing_feature_is_written_as_a_join(self) -> None:
        assert location(
            feature(parts=((3980, 5522),)), 5386
        ) == "join(3981..5386,1..136)"


class TestFeatureTable:
    def test_a_row_for_every_feature(self) -> None:
        features = spread()
        assert len(_rows(features)) == len(features)

    def test_the_lane_is_recorded(self) -> None:
        rows = _rows(spread())
        assert all(isinstance(row[FEATURE_COLUMNS.index("lane")], int) for row in rows)

    def test_a_widened_feature_says_how_wide_it_was_drawn(self) -> None:
        rows = _rows([feature(start=500, end=501)], length=5000)
        assert rows[0][FEATURE_COLUMNS.index("widened")] is True
        assert rows[0][FEATURE_COLUMNS.index("drawn_length")] == 10
        assert rows[0][FEATURE_COLUMNS.index("length")] == 1


class TestWindows:
    def test_a_row_for_every_window(self) -> None:
        profile = measure(sequence(1000), circular=True)
        assert len(window_rows(profile)) == len(profile.windows)

    def test_the_columns_are_what_the_header_says(self) -> None:
        profile = measure(sequence(1000), circular=True)
        assert len(window_rows(profile)[0]) == len(WINDOW_COLUMNS)


class TestKinds:
    def test_every_known_kind_has_a_row(self) -> None:
        from plasmid_map.palette import GROUPS

        assert len(kind_rows(Palette())) == len(GROUPS)

    def test_a_kind_the_program_has_not_heard_of_is_added(self) -> None:
        rows = kind_rows(Palette(), {"something_new": 3})
        assert any(row[0] == "something_new" for row in rows)

    def test_the_columns_are_what_the_header_says(self) -> None:
        assert len(kind_rows(Palette())[0]) == len(KIND_COLUMNS)


class TestSummary:
    def test_every_figure_in_the_summary_is_labelled(self) -> None:
        held = construct(features=spread(), sequence=sequence())
        drawing = draw(
            held,
            select(held),
            Style(gc=True, width=900),
            measure(held.sequence, True),
        )
        table = dict(summary_rows(drawing))
        assert table["length"] == 1000
        assert table["features_drawn"] == drawing.drawn
        assert table["gc_mean"] is not None

    def test_the_counts_add_up_to_what_was_drawn(self) -> None:
        held = construct(features=spread())
        drawing = draw(held, select(held), Style())
        assert sum(counts(drawing).values()) == drawing.drawn
