"""The linear track, including what happens when a circle is cut."""

from __future__ import annotations

import pytest
from conftest import construct, feature, sequence, spread

from plasmid_map.composition import measure
from plasmid_map.drawing import DrawingError, Style, select
from plasmid_map.feature import Span
from plasmid_map.linear import block, draw, pieces_for, scale


def _draw(held, **style):
    return draw(held, select(held), Style(**style))


class TestScale:
    def test_the_first_base_is_at_the_left(self) -> None:
        assert scale(0, 1000, 20.0, 800.0) == 20.0

    def test_the_last_base_is_at_the_right(self) -> None:
        assert scale(1000, 1000, 20.0, 800.0) == 820.0

    def test_the_middle_is_in_the_middle(self) -> None:
        assert scale(500, 1000, 0.0, 800.0) == 400.0

    def test_a_construct_with_no_length_has_no_scale(self) -> None:
        with pytest.raises(DrawingError):
            scale(0, 0, 0.0, 100.0)


class TestCutting:
    def test_a_span_inside_the_construct_is_left_alone(self) -> None:
        assert pieces_for([Span(10, 20)], 1000) == [(Span(10, 20), True)]

    def test_a_span_past_the_end_becomes_the_two_pieces_it_is(self) -> None:
        """In the order the feature reads, tail first, with no line between:
        they are continuous through the cut, not through the middle."""
        assert pieces_for([Span(950, 1050)], 1000) == [
            (Span(950, 1000), False),
            (Span(0, 50), True),
        ]

    def test_the_pieces_hold_the_same_bases(self) -> None:
        pieces = pieces_for([Span(950, 1050)], 1000)
        assert sum(piece.length for piece, _bridge in pieces) == 100

    def test_a_widened_span_past_a_linear_end_is_clamped_not_wrapped(self) -> None:
        """A one base feature at the last base, drawn ten bases wide, is not a
        feature at the start of the construct as well."""
        assert pieces_for([Span(4994, 5004)], 5000, circular=False) == [
            (Span(4994, 5000), True)
        ]

    def test_no_join_line_is_drawn_across_the_cut(self) -> None:
        held = construct(features=(feature(name="wrap", start=950, end=1050),))
        drawing = _draw(held)
        assert drawing.svg.count('class="feature"') == 2
        assert 'class="join"' not in drawing.svg

    def test_a_plus_strand_feature_points_at_the_end_it_really_ends_at(
        self,
    ) -> None:
        """join(951..1000,1..50) ends at base 50, at the left-hand edge."""
        held = construct(
            features=(feature(name="wrap", start=950, end=1050, strand=1),)
        )
        drawing = _draw(held, labels="none")
        import re

        shapes = re.findall(r'<path d="(M[^"]+)"[^>]*class="feature"', drawing.svg)
        assert len(shapes) == 2
        pointed = [shape for shape in shapes if shape.count("L") == 4]
        assert len(pointed) == 1
        # The pointed block is the one at the left, which holds bases 1..50.
        assert float(pointed[0][1:].split()[0]) < drawing.width / 2

    def test_a_crossing_feature_is_drawn_at_both_ends(self) -> None:
        held = construct(features=(feature(name="wrap", start=950, end=1050),))
        drawing = _draw(held)
        assert drawing.svg.count('class="feature"') == 2
        assert len(drawing.split) == 1

    def test_a_terminal_feature_of_a_linear_construct_stays_at_the_end(
        self,
    ) -> None:
        held = construct(
            length=5000,
            topology="linear",
            features=(feature(kind="misc_feature", name="tip", start=4999, end=5000),),
        )
        drawing = _draw(held)
        assert drawing.svg.count('class="feature"') == 1
        assert drawing.split == ()

    def test_the_split_is_reported(self) -> None:
        held = construct(features=(feature(name="wrap", start=950, end=1050),))
        assert [item.name for item in _draw(held).split] == ["wrap"]


class TestBlocks:
    def test_a_plus_strand_block_points_right(self) -> None:
        drawn = block(0.0, 100.0, 0.0, 15.0, 1, "#000000", "x")
        backwards = block(0.0, 100.0, 0.0, 15.0, -1, "#000000", "x")
        assert drawn != backwards

    def test_a_block_with_no_strand_is_a_rectangle(self) -> None:
        drawn = block(0.0, 100.0, 0.0, 15.0, 0, "#000000", "x")
        assert drawn.count("L") == 3

    def test_an_arrowhead_never_runs_past_the_block(self) -> None:
        drawn = block(0.0, 2.0, 0.0, 15.0, 1, "#000000", "x")
        assert "-" not in drawn.split('d="')[1].split('"')[0]

    def test_a_block_carries_its_position_as_a_tooltip(self) -> None:
        assert "<title>" in block(0.0, 10.0, 0.0, 15.0, 1, "#000000", "bla 1..10")


class TestContents:
    def test_every_feature_is_drawn(self) -> None:
        held = construct(topology="linear", features=spread()[:-1])
        drawing = _draw(held)
        assert drawing.drawn == len(drawing.selection.kept)

    def test_every_feature_has_a_label(self) -> None:
        held = construct(topology="linear", features=spread()[:-1])
        drawing = _draw(held)
        assert len(drawing.placed.labels) == drawing.drawn

    def test_no_two_labels_overlap(self) -> None:
        held = construct(topology="linear", features=spread()[:-1])
        assert _draw(held).label_collisions == []

    def test_the_axis_is_drawn(self) -> None:
        assert 'class="axis"' in _draw(construct(features=spread())).svg

    def test_the_heading_carries_the_name_and_the_topology(self) -> None:
        held = construct(topology="linear", name="cassette", features=spread()[:-1])
        drawing = _draw(held)
        assert "cassette" in drawing.svg
        assert "1,000 bp, linear" in drawing.svg

    def test_spliced_parts_are_joined_by_a_line(self) -> None:
        held = construct(
            topology="linear", features=(feature(parts=((10, 100), (300, 400))),)
        )
        assert 'class="join"' in _draw(held).svg

    def test_labels_can_be_left_off(self) -> None:
        drawing = _draw(construct(features=spread()), labels="none")
        assert drawing.placed is None


class TestTracks:
    def test_a_gc_track_is_drawn_below_the_ruler(self) -> None:
        held = construct(
            topology="linear", features=spread()[:-1], sequence=sequence()
        )
        drawing = draw(
            held,
            select(held),
            Style(gc=True),
            measure(held.sequence, False),
        )
        assert 'class="gc-fill"' in drawing.svg

    def test_a_gc_track_needs_the_sequence(self) -> None:
        held = construct(topology="linear", features=spread()[:-1])
        with pytest.raises(DrawingError):
            draw(held, select(held), Style(gc=True), None)

    def test_a_skew_track_is_drawn_when_one_is_asked_for(self) -> None:
        """Asking for a track and getting a map without one is the worst kind
        of success."""
        held = construct(
            topology="linear", features=spread()[:-1], sequence=sequence()
        )
        drawing = draw(
            held, select(held), Style(skew=True), measure(held.sequence, False)
        )
        assert drawing.svg.count('class="skew-up"') > 0
        assert drawing.svg.count('class="skew-down"') > 0

    def test_a_skew_track_alone_needs_the_sequence(self) -> None:
        held = construct(topology="linear", features=spread()[:-1])
        with pytest.raises(DrawingError):
            draw(held, select(held), Style(skew=True), None)

    def test_both_tracks_are_drawn_when_both_are_asked_for(self) -> None:
        held = construct(
            topology="linear", features=spread()[:-1], sequence=sequence()
        )
        drawing = draw(
            held,
            select(held),
            Style(gc=True, skew=True),
            measure(held.sequence, False),
        )
        assert 'class="gc-fill"' in drawing.svg
        assert 'class="skew-up"' in drawing.svg


class TestRefusals:
    def test_rotating_a_linear_track_is_refused_rather_than_ignored(self) -> None:
        held = construct(topology="linear", features=spread()[:-1])
        with pytest.raises(DrawingError) as raised:
            draw(held, select(held), Style(start=100))
        assert "circular" in str(raised.value)

    def test_a_track_with_no_room_to_draw_in_is_refused(self) -> None:
        held = construct(topology="linear", features=spread()[:-1])
        with pytest.raises(DrawingError):
            draw(held, select(held), Style(width=50))


class TestDeterminism:
    def test_the_same_construct_draws_the_same_bytes(self) -> None:
        held = construct(
            topology="linear", features=spread()[:-1], sequence=sequence()
        )
        profile = measure(held.sequence, False)
        style = Style(gc=True, legend=True)
        assert (
            draw(held, select(held), style, profile).svg
            == draw(held, select(held), style, profile).svg
        )


class TestRealRecords:
    def test_aav2_draws_its_spliced_genes_with_a_join_line(self, aav2) -> None:
        selection = select(aav2, kinds=["CDS"])
        drawing = draw(aav2, selection, Style(width=900))
        assert 'class="join"' in drawing.svg
        assert drawing.split == ()

    def test_a_circular_record_drawn_as_a_line_is_cut_at_the_origin(
        self, phix174
    ) -> None:
        selection = select(phix174, kinds=["gene"])
        drawing = draw(phix174, selection, Style(width=900))
        assert len(drawing.split) == 3
        assert drawing.label_collisions == []
