"""Lane assignment: nothing in one lane may overlap anything else in it."""

from __future__ import annotations

import pytest
from conftest import feature, spread

from plasmid_map.feature import Span, minimum_width, overlapping
from plasmid_map.lanes import assign, deepest, gap_for, lane_radii, lane_rows


def _lanes_are_clear(layout) -> bool:
    """Whether every lane holds features that keep out of each other's way."""
    for lane in range(layout.lanes):
        held = layout.in_lane(lane)
        for first in range(len(held)):
            for second in range(first + 1, len(held)):
                if overlapping(
                    held[first].extent,
                    held[second].extent,
                    layout.length,
                    layout.circular,
                ):
                    return False
    return True


class TestAssignment:
    def test_nothing_in_a_lane_overlaps_anything_else_in_it(self) -> None:
        layout = assign(spread(), 1000, True, minimum_width(1000))
        assert _lanes_are_clear(layout)

    def test_every_feature_is_placed(self) -> None:
        features = spread()
        layout = assign(features, 1000, True, minimum_width(1000))
        assert len(layout.placements) == len(features)

    def test_an_enclosing_feature_takes_the_outer_lane(self) -> None:
        """A gene outside the coding sequence it holds is what a reader wants."""
        layout = assign(
            [
                feature(kind="gene", name="gene", start=0, end=900),
                feature(kind="CDS", name="cds", start=10, end=800),
            ],
            1000,
            True,
            1,
        )
        placed = {item.feature.name: item.lane for item in layout.placements}
        assert placed["gene"] == 0
        assert placed["cds"] == 1

    def test_features_that_do_not_meet_share_a_lane(self) -> None:
        layout = assign(
            [feature(start=0, end=100), feature(start=400, end=500)], 1000, True, 1
        )
        assert layout.lanes == 1

    def test_a_feature_near_the_end_and_one_near_the_start_are_neighbours(
        self,
    ) -> None:
        """Only when the construct is circular."""
        near = [feature(name="a", start=950, end=1000), feature(name="b", start=0, end=40)]
        assert assign(near, 1000, True, 1, gap=20).lanes == 2
        assert assign(near, 1000, False, 1, gap=20).lanes == 1

    def test_a_feature_crossing_the_origin_blocks_its_own_lane_at_both_ends(
        self,
    ) -> None:
        layout = assign(
            [
                feature(name="wrap", start=980, end=1030),
                feature(name="early", start=10, end=20),
            ],
            1000,
            True,
            1,
        )
        assert layout.lanes == 2

    def test_lanes_are_counted_as_they_are_used(self) -> None:
        layout = assign(spread(), 1000, True, minimum_width(1000))
        assert layout.lanes == deepest(layout.placements)


class TestWidening:
    def test_a_one_base_feature_is_drawn_wider_than_it_is(self) -> None:
        layout = assign([feature(start=500, end=501)], 5000, True, 10)
        assert layout.placements[0].spans[0].length == 10
        assert layout.placements[0].widened

    def test_a_long_feature_is_not_widened(self) -> None:
        layout = assign([feature(start=0, end=900)], 5000, True, 10)
        assert not layout.placements[0].widened

    def test_widened_features_are_kept_apart_as_drawn_not_as_measured(self) -> None:
        """Two marks three bases apart are drawn on top of each other unless
        the lanes know how wide they were drawn."""
        marks = [
            feature(name="a", start=500, end=501),
            feature(name="b", start=503, end=504),
        ]
        layout = assign(marks, 5000, True, 20)
        assert layout.lanes == 2

    def test_which_ones_were_widened_is_reported(self) -> None:
        layout = assign(
            [feature(name="tiny", start=10, end=11), feature(name="big", start=100, end=900)],
            5000,
            True,
            10,
        )
        assert [item.feature.name for item in layout.widened] == ["tiny"]


class TestExtent:
    def test_the_extent_covers_every_part(self) -> None:
        layout = assign([feature(parts=((10, 50), (200, 260)))], 1000, True, 1)
        assert layout.placements[0].extent == Span(10, 260)


class TestGeometry:
    def test_lane_zero_is_the_outermost_ring(self) -> None:
        radii = lane_radii(200.0, 15.0, 4.0, 3)
        assert radii[0][1] == 200.0
        assert radii[1][1] < radii[0][1]

    def test_rings_do_not_touch(self) -> None:
        radii = lane_radii(200.0, 15.0, 4.0, 3)
        for outer, inner in zip(radii, radii[1:]):
            assert inner[1] < outer[0]

    def test_the_first_row_is_at_the_top(self) -> None:
        rows = lane_rows(40.0, 15.0, 4.0, 3)
        assert rows[0][0] == 40.0
        assert rows[1][0] > rows[0][0]

    def test_one_lane_still_gets_a_band(self) -> None:
        assert len(lane_radii(200.0, 15.0, 4.0, 0)) == 1

    def test_the_gap_is_at_least_one_base(self) -> None:
        assert gap_for(10) == 1


class TestRealRecords:
    @pytest.mark.parametrize("name", ["pbr322", "phix174", "aav2"])
    def test_a_real_record_assigns_into_clear_lanes(self, name, request) -> None:
        construct = request.getfixturevalue(name)
        layout = assign(
            construct.features,
            construct.length,
            construct.is_circular,
            minimum_width(construct.length),
        )
        assert _lanes_are_clear(layout)
        assert len(layout.placements) == len(construct.features)
