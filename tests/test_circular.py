"""The circular map: what it holds, what it reports, and what it refuses."""

from __future__ import annotations

import re

import pytest
from conftest import construct, feature, sequence, spread

from plasmid_map.circular import (
    BAND_SHARE,
    LEAST_CENTRE,
    centre_reserve,
    choose_ticks,
    draw,
    fit_bands,
    side_room,
)
from plasmid_map.composition import measure
from plasmid_map.drawing import DrawingError, Style, select
from plasmid_map.geometry import Point
from plasmid_map.labels import Request


def _draw(held, **style):
    return draw(held, select(held), Style(**style))


class TestContents:
    def test_every_selected_feature_is_drawn(self) -> None:
        held = construct(features=spread())
        drawing = _draw(held)
        assert drawing.drawn == len(drawing.selection.kept)

    def test_every_feature_has_a_shape_in_the_file(self) -> None:
        held = construct(features=spread())
        drawing = _draw(held)
        assert drawing.svg.count('class="feature"') >= drawing.drawn

    def test_every_feature_has_a_label(self) -> None:
        held = construct(features=spread())
        drawing = _draw(held)
        assert len(drawing.placed.labels) == drawing.drawn

    def test_no_two_labels_overlap(self) -> None:
        held = construct(features=spread())
        assert _draw(held).label_collisions == []

    def test_labels_can_be_left_off(self) -> None:
        held = construct(features=spread())
        drawing = _draw(held, labels="none")
        assert drawing.placed is None
        assert 'class="label"' not in drawing.svg

    def test_the_backbone_is_drawn(self) -> None:
        assert 'class="backbone"' in _draw(construct(features=spread())).svg

    def test_the_name_and_the_length_are_in_the_middle(self) -> None:
        held = construct(features=spread(), name="pMAP1")
        drawing = _draw(held)
        assert "pMAP1" in drawing.svg
        assert "1,000 bp" in drawing.svg

    def test_a_title_replaces_the_name(self) -> None:
        held = construct(features=spread(), name="pMAP1")
        assert "given name" in _draw(held, title="given name").svg

    def test_a_linear_record_drawn_as_a_circle_says_so(self) -> None:
        held = construct(topology="linear", features=(feature(start=0, end=100),))
        assert "the record calls it linear" in _draw(held).svg


class TestOrigin:
    def test_a_feature_across_the_origin_is_one_arc(self) -> None:
        held = construct(features=(feature(name="wrap", start=950, end=1050),))
        drawing = _draw(held)
        assert len(drawing.crossing) == 1
        shapes = re.findall(r'class="feature"', drawing.svg)
        assert len(shapes) == 1

    def test_its_arc_is_not_drawn_the_long_way_round(self) -> None:
        held = construct(features=(feature(name="wrap", start=950, end=1050),))
        drawing = _draw(held, labels="none")
        found = re.search(r'd="(M[^"]+)"[^>]*class="feature"', drawing.svg)
        assert found is not None
        # A 100 base span of a 1,000 base plasmid is 36 degrees, so the arc
        # must not carry the large-arc flag.
        assert " 1 1 " not in found.group(1)


class TestBands:
    def test_bands_are_left_alone_when_they_fit(self) -> None:
        style = Style()
        assert fit_bands(400.0, 2, style) == (style.track, style.separation)

    def test_bands_are_thinned_when_there_are_too_many_lanes(self) -> None:
        style = Style()
        track, separation = fit_bands(100.0, 12, style)
        assert track < style.track
        assert separation < style.separation

    def test_thinned_bands_keep_to_their_share_of_the_radius(self) -> None:
        style = Style()
        radius = 120.0
        lanes = 10
        track, separation = fit_bands(radius, lanes, style)
        needed = (lanes - 1) * (track + separation) + track
        assert needed <= radius * BAND_SHARE + 0.01

    def test_the_thinning_is_reported(self) -> None:
        many = tuple(
            feature(name=f"f{number}", start=0, end=900)
            for number in range(12)
        )
        drawing = _draw(construct(features=many), width=600)
        assert drawing.thinned
        assert drawing.track < drawing.style.track


class TestCentre:
    def test_a_long_name_reserves_more_of_the_middle(self) -> None:
        style = Style()
        assert centre_reserve("a very long construct name indeed", "1 bp", style) > (
            centre_reserve("p1", "1 bp", style)
        )

    def test_there_is_always_some_room_kept(self) -> None:
        assert centre_reserve("p", "1 bp", Style()) == LEAST_CENTRE

    def test_a_track_that_would_reach_the_name_is_refused(self) -> None:
        held = construct(features=spread(), sequence=sequence())
        profile = measure(held.sequence, True)
        with pytest.raises(DrawingError) as raised:
            draw(held, select(held), Style(width=420, gc=True, skew=True), profile)
        assert "kept clear in the middle" in str(raised.value)


class TestTracks:
    def test_a_gc_track_needs_the_sequence(self) -> None:
        held = construct(features=spread())
        with pytest.raises(DrawingError) as raised:
            draw(held, select(held), Style(gc=True), None)
        assert "carries none" in str(raised.value)

    def test_a_gc_track_is_drawn_when_there_is_one(self) -> None:
        held = construct(features=spread(), sequence=sequence())
        drawing = draw(
            held, select(held), Style(gc=True, width=900), measure(held.sequence, True)
        )
        assert 'class="gc-fill"' in drawing.svg
        assert drawing.profile is not None

    def test_a_skew_track_is_two_coloured(self) -> None:
        held = construct(features=spread(), sequence=sequence())
        drawing = draw(
            held,
            select(held),
            Style(skew=True, width=900),
            measure(held.sequence, True),
        )
        assert 'class="skew-up"' in drawing.svg
        assert 'class="skew-down"' in drawing.svg

    def test_the_gc_fill_lies_between_the_curve_and_the_mean(self) -> None:
        """Which is the even-odd rule, and has to be set or the ring fills."""
        held = construct(features=spread(), sequence=sequence())
        drawing = draw(
            held, select(held), Style(gc=True, width=900), measure(held.sequence, True)
        )
        assert "fill-rule: evenodd" in drawing.svg


class TestRuler:
    def test_the_ruler_can_be_left_off(self) -> None:
        drawing = _draw(construct(features=spread()), ruler=False)
        assert 'class="tick-label"' not in drawing.svg

    def test_the_tick_labels_do_not_collide(self) -> None:
        marks, _step = choose_ticks(
            construct(features=()), Point(0.0, 0.0), 60.0, Style()
        )
        from plasmid_map.circular import _ticks_fit

        assert _ticks_fit(marks, Point(0.0, 0.0), 60.0, Style())

    def test_a_small_map_is_given_a_coarser_ruler(self) -> None:
        held = construct(length=100_000, features=(feature(start=0, end=100),))
        _fine, coarse = choose_ticks(held, Point(0.0, 0.0), 40.0, Style())
        _wide, fine = choose_ticks(held, Point(0.0, 0.0), 400.0, Style())
        assert coarse >= fine

    def test_the_step_used_is_reported(self) -> None:
        assert _draw(construct(features=spread())).tick_step > 0


class TestCanvas:
    def test_the_canvas_is_the_width_that_was_asked_for(self) -> None:
        drawing = _draw(construct(features=spread()), width=820)
        assert drawing.width == pytest.approx(820.0)

    def test_the_canvas_grows_tall_enough_to_hold_every_label(self) -> None:
        crowded = tuple(
            feature(name=f"feature number {number}", start=number, end=number + 2)
            for number in range(40)
        )
        drawing = _draw(construct(features=crowded))
        bounds = drawing.placed.bounds()
        assert bounds.height <= drawing.height

    def test_a_map_with_no_room_for_its_labels_is_refused_with_advice(self) -> None:
        held = construct(features=spread())
        with pytest.raises(DrawingError) as raised:
            _draw(held, width=120)
        assert "--width" in str(raised.value)


class TestDeterminism:
    def test_the_same_construct_draws_the_same_bytes(self) -> None:
        held = construct(features=spread(), sequence=sequence())
        profile = measure(held.sequence, True)
        style = Style(gc=True, width=900, legend=True)
        first = draw(held, select(held), style, profile).svg
        second = draw(held, select(held), style, profile).svg
        assert first == second

    def test_nothing_in_the_file_depends_on_the_clock(self) -> None:
        drawing = _draw(construct(features=spread()))
        assert not re.search(r"20\d\d-\d\d-\d\d", drawing.svg)


class TestLegend:
    def test_a_legend_names_the_groups_present(self) -> None:
        held = construct(
            features=(
                feature(kind="CDS", start=0, end=100),
                feature(kind="rep_origin", start=200, end=300),
            )
        )
        drawing = _draw(held, legend=True)
        assert "coding" in drawing.svg
        assert "origin" in drawing.svg

    def test_there_is_no_legend_unless_one_is_asked_for(self) -> None:
        drawing = _draw(construct(features=spread()))
        assert 'class="legend-swatch"' not in drawing.svg


class TestSideRoom:
    def test_a_label_on_the_right_asks_for_room_on_the_right(self) -> None:
        left, right = side_room([Request("a long label", 90.0)], Style())
        assert right > 0 and left == 0

    def test_a_label_on_the_left_asks_for_room_on_the_left(self) -> None:
        left, right = side_room([Request("a long label", 270.0)], Style())
        assert left > 0 and right == 0


class TestRealRecords:
    def test_pbr322_draws_with_nothing_lost(self, pbr322) -> None:
        selection = select(pbr322)
        drawing = draw(pbr322, selection, Style(width=900))
        assert drawing.drawn == len(selection.kept)
        assert drawing.label_collisions == []

    def test_phix174_draws_its_three_crossing_genes_as_single_arcs(
        self, phix174
    ) -> None:
        selection = select(phix174, kinds=["CDS", "gene"])
        drawing = draw(phix174, selection, Style(width=900))
        assert len(drawing.crossing) == 6
        for item in drawing.layout.placements:
            assert len(item.spans) == 1

    def test_the_source_feature_is_left_out_by_default(self, pbr322) -> None:
        drawing = draw(pbr322, select(pbr322), Style(width=900))
        assert drawing.selection.whole_record
        assert all(
            item.feature.kind != "source" for item in drawing.layout.placements
        )
