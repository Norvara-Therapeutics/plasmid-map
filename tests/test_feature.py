"""The coordinate model, which every other result depends on."""

from __future__ import annotations

import pytest
from conftest import construct, feature

from plasmid_map.feature import (
    Construct,
    Feature,
    FeatureError,
    Span,
    drawn_spans,
    joined,
    minimum_width,
    overlapping,
    sorted_features,
    strand_of,
    within,
)


class TestSpan:
    def test_a_span_is_half_open(self) -> None:
        assert Span(85, 1276).length == 1191

    def test_a_site_between_two_bases_has_no_length(self) -> None:
        """GenBank writes those ``102^103``; they are marks, not stretches."""
        assert Span(102, 102).length == 0

    def test_a_span_cannot_start_before_the_first_base(self) -> None:
        with pytest.raises(FeatureError):
            Span(-1, 10)

    def test_a_span_cannot_end_before_it_starts(self) -> None:
        with pytest.raises(FeatureError):
            Span(10, 9)

    def test_an_end_past_the_construct_is_allowed(self) -> None:
        """It is how a feature crossing the origin is held."""
        assert Span(3980, 5522).length == 1542

    def test_widening_grows_a_span_about_its_middle(self) -> None:
        widened = Span(500, 501).widened(21)
        assert widened.length == 21
        assert widened.start <= 500 < widened.end

    def test_widening_leaves_a_long_span_alone(self) -> None:
        assert Span(0, 100).widened(21) == Span(0, 100)

    def test_widening_never_runs_off_the_front(self) -> None:
        assert Span(0, 1).widened(21).start == 0

    def test_shifting_moves_both_ends(self) -> None:
        assert Span(10, 20).shifted(5) == Span(15, 25)


class TestOverlap:
    def test_spans_that_touch_do_not_overlap(self) -> None:
        assert not overlapping(Span(0, 10), Span(10, 20), 1000, True)

    def test_spans_that_share_a_base_overlap(self) -> None:
        assert overlapping(Span(0, 11), Span(10, 20), 1000, True)

    def test_a_span_near_the_end_meets_one_near_the_start(self) -> None:
        """Only on a circular construct; the two are neighbours there."""
        assert overlapping(Span(990, 1000), Span(0, 5), 1000, True, padding=5)
        assert not overlapping(Span(990, 1000), Span(0, 5), 1000, False, padding=5)

    def test_the_turn_is_taken_in_arithmetic_not_by_building_a_span(self) -> None:
        """A span shifted back a turn starts before base 1, which is no place."""
        assert overlapping(Span(0, 10), Span(995, 1005), 1000, True)

    def test_padding_keeps_neighbours_apart(self) -> None:
        assert not overlapping(Span(0, 10), Span(12, 20), 1000, True)
        assert overlapping(Span(0, 10), Span(12, 20), 1000, True, padding=5)


class TestFeature:
    def test_its_length_is_the_sum_of_its_parts(self) -> None:
        """Not the distance from its first coordinate to its last."""
        spliced = feature(parts=((320, 1906), (2227, 2252)))
        assert spliced.length == 1586 + 25
        assert spliced.end - spliced.start == 1932

    def test_an_origin_crossing_feature_has_the_length_of_its_bases(self) -> None:
        crossing = feature(parts=((3980, 5522),))
        assert crossing.length == 1542
        assert crossing.crosses_origin(5386)

    def test_parts_are_put_in_order(self) -> None:
        out_of_order = Feature(
            kind="CDS", parts=(Span(500, 600), Span(100, 200)), name="x"
        )
        assert [part.start for part in out_of_order.parts] == [100, 500]

    def test_a_feature_needs_a_kind(self) -> None:
        with pytest.raises(FeatureError):
            Feature(kind="", parts=(Span(0, 10),))

    def test_a_feature_needs_a_position(self) -> None:
        with pytest.raises(FeatureError):
            Feature(kind="CDS", parts=())

    def test_a_strand_is_a_sign(self) -> None:
        with pytest.raises(FeatureError):
            Feature(kind="CDS", parts=(Span(0, 10),), strand=2)

    def test_the_label_falls_back_to_the_kind(self) -> None:
        assert Feature(kind="rep_origin", parts=(Span(0, 1),)).label == "rep_origin"

    def test_a_qualifier_that_is_not_there_is_none(self) -> None:
        held = Feature(
            kind="CDS", parts=(Span(0, 9),), qualifiers=(("gene", "bla"),)
        )
        assert held.qualifier("gene") == "bla"
        assert held.qualifier("product") is None

    def test_partial_ends_are_kept(self) -> None:
        marked = Feature(kind="CDS", parts=(Span(0, 10),), partial=(True, False))
        assert marked.is_partial

    def test_with_parts_keeps_everything_else(self) -> None:
        original = Feature(
            kind="CDS",
            name="bla",
            parts=(Span(0, 10),),
            strand=-1,
            partial=(True, False),
            qualifiers=(("gene", "bla"),),
        )
        moved = original.with_parts([Span(20, 30)])
        assert moved.strand == -1
        assert moved.partial == (True, False)
        assert moved.qualifiers == (("gene", "bla"),)


class TestJoined:
    def test_parts_meeting_at_the_origin_become_one(self) -> None:
        parts = joined((Span(3980, 5386), Span(0, 136)), 5386, True)
        assert parts == (Span(3980, 5522),)

    def test_a_join_with_a_gap_stays_two_parts(self) -> None:
        """``join(4000..5386,50..100)`` does not meet at the origin."""
        parts = joined((Span(3999, 5386), Span(49, 100)), 5386, True)
        assert len(parts) == 2

    def test_a_linear_molecule_never_joins_at_the_origin(self) -> None:
        parts = joined((Span(3980, 5386), Span(0, 136)), 5386, False)
        assert len(parts) == 2

    def test_a_spliced_feature_is_left_alone(self) -> None:
        parts = joined((Span(320, 1906), Span(2227, 2252)), 4679, False)
        assert len(parts) == 2


class TestWithin:
    def test_a_feature_past_a_linear_end_is_not_within_it(self) -> None:
        assert not within(feature(start=900, end=1100), 1000, False)

    def test_a_feature_past_a_circular_end_is_within_it(self) -> None:
        assert within(feature(start=900, end=1100), 1000, True)

    def test_a_feature_may_not_go_more_than_once_round(self) -> None:
        assert not within(feature(start=10, end=1200), 1000, True)


class TestConstruct:
    def test_a_sequence_has_to_be_as_long_as_the_record_says(self) -> None:
        with pytest.raises(FeatureError) as raised:
            Construct(name="x", length=100, sequence="ACGT")
        assert "100" in str(raised.value) and "4" in str(raised.value)

    def test_a_feature_off_the_end_is_refused_with_its_position(self) -> None:
        with pytest.raises(FeatureError) as raised:
            construct(
                length=100, topology="linear", features=(feature(start=90, end=200),)
            )
        assert "91..200" in str(raised.value)

    def test_not_a_topology_is_refused_by_name(self) -> None:
        with pytest.raises(FeatureError) as raised:
            Construct(name="x", length=10, topology="round")
        assert "circular" in str(raised.value)

    def test_kinds_are_counted_most_common_first(self) -> None:
        held = construct(
            features=(
                feature(kind="CDS", start=0, end=10),
                feature(kind="CDS", start=20, end=30),
                feature(kind="gene", start=40, end=50),
            )
        )
        assert held.kinds() == [("CDS", 2), ("gene", 1)]

    def test_crossing_finds_the_features_past_the_origin(self) -> None:
        held = construct(
            length=1000, features=(feature(start=980, end=1020), feature(start=0, end=5))
        )
        assert len(held.crossing()) == 1


class TestOrdering:
    def test_the_longest_of_two_starting_together_comes_first(self) -> None:
        """So an enclosing feature takes the outer lane."""
        ordered = sorted_features(
            [feature(name="short", start=0, end=10), feature(name="long", start=0, end=99)]
        )
        assert [item.name for item in ordered] == ["long", "short"]


class TestMinimumWidth:
    def test_a_short_construct_still_gets_a_whole_base(self) -> None:
        assert minimum_width(100) == 1

    def test_it_is_a_small_share_of_a_long_one(self) -> None:
        assert minimum_width(5000) == 10

    def test_drawn_spans_widen_every_part(self) -> None:
        spans = drawn_spans(feature(parts=((10, 11), (500, 501))), 5000, 10)
        assert [span.length for span in spans] == [10, 10]

    def test_drawn_spans_leave_a_long_part_alone(self) -> None:
        spans = drawn_spans(feature(start=0, end=900), 5000, 10)
        assert spans[0].length == 900


class TestStrands:
    @pytest.mark.parametrize(
        "written,expected",
        [("+", 1), ("-", -1), (".", 0), ("", 0), ("?", 0), ("1", 1), ("-1", -1)],
    )
    def test_a_strand_is_read_from_how_a_format_writes_it(
        self, written: str, expected: int
    ) -> None:
        assert strand_of(written) == expected

    def test_anything_else_is_refused(self) -> None:
        with pytest.raises(FeatureError):
            strand_of("forwards")
