"""Angles, arcs and the numbers they are written with."""

from __future__ import annotations

import math

import pytest

from plasmid_map.geometry import (
    Box,
    GeometryError,
    Point,
    angle_at,
    arc_band,
    arc_line,
    chord_length,
    number,
    polar,
    polyline,
    radial,
    ring,
    unwrapped_angles,
)

CENTRE = Point(0.0, 0.0)


class TestNumbers:
    def test_a_coordinate_is_written_to_two_places(self) -> None:
        assert number(1 / 3) == "0.33"

    def test_negative_zero_is_written_as_zero(self) -> None:
        """Otherwise the same map could produce two different files."""
        assert number(-0.0) == "0.00"
        assert number(-0.001) == "0.00"

    def test_a_whole_number_keeps_its_places(self) -> None:
        assert number(12) == "12.00"


class TestAngles:
    def test_the_origin_is_at_twelve_oclock(self) -> None:
        assert angle_at(0, 1000) == 0.0

    def test_angles_increase_clockwise(self) -> None:
        assert angle_at(250, 1000) == 90.0

    def test_start_turns_the_map(self) -> None:
        assert angle_at(250, 1000, start=250) == 0.0

    def test_a_position_past_the_end_comes_round(self) -> None:
        assert angle_at(1250, 1000) == 90.0

    def test_a_construct_with_no_length_has_no_angles(self) -> None:
        with pytest.raises(GeometryError):
            angle_at(0, 0)

    def test_a_span_crossing_twelve_oclock_is_not_drawn_the_long_way(self) -> None:
        first, last = unwrapped_angles(900, 1100, 1000)
        assert first == pytest.approx(324.0)
        assert last == pytest.approx(396.0)
        assert last > first


class TestPolar:
    def test_twelve_oclock_is_straight_up(self) -> None:
        point = polar(CENTRE, 100.0, 0.0)
        assert point.x == pytest.approx(0.0)
        assert point.y == pytest.approx(-100.0)

    def test_ninety_degrees_is_to_the_right(self) -> None:
        point = polar(CENTRE, 100.0, 90.0)
        assert point.x == pytest.approx(100.0)
        assert point.y == pytest.approx(0.0, abs=1e-9)

    def test_every_point_is_on_the_circle(self) -> None:
        for angle in range(0, 360, 17):
            point = polar(CENTRE, 42.0, float(angle))
            assert math.hypot(point.x, point.y) == pytest.approx(42.0)

    def test_a_point_writes_as_a_pair(self) -> None:
        assert str(Point(1.0, -2.0)) == "1.00,-2.00"


class TestBands:
    def test_a_band_needs_an_outer_radius_larger_than_its_inner(self) -> None:
        with pytest.raises(GeometryError):
            arc_band(CENTRE, 100.0, 100.0, 0.0, 90.0)

    def test_a_band_is_a_closed_path(self) -> None:
        assert arc_band(CENTRE, 90.0, 100.0, 0.0, 90.0).endswith("Z")

    def test_a_band_of_more_than_half_a_turn_sets_the_large_arc_flag(self) -> None:
        assert " 1 1 " in arc_band(CENTRE, 90.0, 100.0, 0.0, 200.0)

    def test_a_whole_turn_is_drawn_as_two_rings(self) -> None:
        drawn = arc_band(CENTRE, 90.0, 100.0, 0.0, 360.0)
        assert drawn.count("M") == 2

    def test_an_arrowhead_makes_a_point(self) -> None:
        plain = arc_band(CENTRE, 90.0, 100.0, 0.0, 90.0)
        pointed = arc_band(CENTRE, 90.0, 100.0, 0.0, 90.0, arrow="end", arrow_degrees=10)
        assert plain != pointed

    def test_a_feature_shorter_than_its_arrowhead_becomes_a_triangle(self) -> None:
        """Rather than an arrowhead pointing back past its own start."""
        drawn = arc_band(
            CENTRE, 90.0, 100.0, 0.0, 2.0, arrow="end", arrow_degrees=10.0
        )
        assert "A" not in drawn
        assert drawn.count("L") == 2

    def test_a_backwards_band_is_drawn_the_same_as_a_forwards_one(self) -> None:
        assert arc_band(CENTRE, 90.0, 100.0, 90.0, 0.0) == arc_band(
            CENTRE, 90.0, 100.0, 0.0, 90.0
        )

    def test_the_two_arrow_directions_differ(self) -> None:
        forwards = arc_band(CENTRE, 90.0, 100.0, 0.0, 90.0, "end", 10.0)
        backwards = arc_band(CENTRE, 90.0, 100.0, 0.0, 90.0, "start", 10.0)
        assert forwards != backwards


class TestLines:
    def test_a_ring_closes(self) -> None:
        assert ring(CENTRE, 50.0).endswith("Z")

    def test_an_arc_line_starts_where_it_is_told(self) -> None:
        assert arc_line(CENTRE, 50.0, 0.0, 90.0).startswith("M0.00,-50.00")

    def test_a_whole_turn_of_arc_line_is_a_ring(self) -> None:
        assert arc_line(CENTRE, 50.0, 0.0, 360.0) == ring(CENTRE, 50.0)

    def test_a_radial_line_runs_out_from_the_centre(self) -> None:
        assert radial(CENTRE, 40.0, 50.0, 0.0) == "M0.00,-40.00L0.00,-50.00"

    def test_an_empty_polyline_draws_nothing(self) -> None:
        assert polyline([]) == ""

    def test_a_polyline_visits_every_point(self) -> None:
        drawn = polyline([Point(0, 0), Point(1, 1), Point(2, 2)])
        assert drawn.count("L") == 2


class TestChords:
    def test_half_a_turn_spans_the_diameter(self) -> None:
        assert chord_length(50.0, 180.0) == pytest.approx(100.0)

    def test_a_small_angle_is_near_the_arc_length(self) -> None:
        assert chord_length(100.0, 1.0) == pytest.approx(math.radians(1.0) * 100, rel=1e-4)


class TestBoxes:
    def test_a_box_knows_its_far_edges(self) -> None:
        box = Box(10.0, 20.0, 5.0, 6.0)
        assert box.right == 15.0
        assert box.bottom == 26.0

    def test_growing_adds_a_margin_on_every_side(self) -> None:
        grown = Box(10.0, 10.0, 10.0, 10.0).grown(2.0)
        assert (grown.x, grown.y, grown.width, grown.height) == (8.0, 8.0, 14.0, 14.0)

    def test_boxes_apart_do_not_overlap(self) -> None:
        assert not Box(0, 0, 10, 10).overlaps(Box(20, 20, 10, 10))
