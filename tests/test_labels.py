"""Label placement, which is the thing this program is for.

The claims tested here are the ones the README makes: no label is dropped, no
two labels overlap, labels keep the order their features are in, and each is
moved as little as the arrangement allows.
"""

from __future__ import annotations

import pytest

from plasmid_map.geometry import Box, Point
from plasmid_map.labels import (
    INK_HEIGHT,
    LINE_HEIGHT,
    Along,
    LabelError,
    Request,
    estimate_width,
    place_circular,
    separate,
    shorten,
    stack_linear,
)

CENTRE = Point(0.0, 0.0)


class TestShorten:
    def test_a_short_label_is_left_alone(self) -> None:
        assert shorten("bla", 28) == "bla"

    def test_a_long_label_is_cut_with_an_ellipsis(self) -> None:
        cut = shorten("gamma-delta insertion target sequence", 20)
        assert len(cut) == 20
        assert cut.endswith("\u2026")

    def test_whitespace_is_collapsed_so_a_label_is_one_line(self) -> None:
        assert shorten("one\ntwo   three") == "one two three"


class TestSeparate:
    def test_nothing_moves_when_nothing_collides(self) -> None:
        assert separate([0.0, 50.0, 100.0], 10.0) == [0.0, 50.0, 100.0]

    def test_everything_ends_up_at_least_a_spacing_apart(self) -> None:
        settled = separate([0.0, 0.0, 0.0, 0.0], 10.0)
        gaps = [second - first for first, second in zip(settled, settled[1:])]
        assert all(gap >= 10.0 - 1e-9 for gap in gaps)

    def test_a_block_is_centred_on_what_its_members_wanted(self) -> None:
        settled = separate([0.0, 0.0, 0.0], 10.0)
        assert sum(settled) == pytest.approx(0.0)

    def test_the_order_is_kept(self) -> None:
        """A label whose feature is higher stays higher, so a leader works."""
        wanted = [30.0, 10.0, 20.0, 11.0]
        settled = separate(wanted, 8.0)
        by_wish = sorted(range(len(wanted)), key=lambda index: wanted[index])
        places = [settled[index] for index in by_wish]
        assert places == sorted(places)

    def test_nothing_is_dropped(self) -> None:
        assert len(separate([1.0] * 40, 5.0)) == 40

    def test_an_empty_list_settles_to_an_empty_list(self) -> None:
        assert separate([], 10.0) == []

    def test_labels_need_room_between_them(self) -> None:
        with pytest.raises(LabelError):
            separate([0.0, 1.0], 0.0)

    def test_a_long_run_spreads_symmetrically_about_its_wish(self) -> None:
        settled = separate([0.0] * 5, 10.0)
        assert settled[0] == pytest.approx(-20.0)
        assert settled[-1] == pytest.approx(20.0)


class TestCircular:
    def _requests(self, count: int) -> list[Request]:
        return [
            Request(text=f"feature {number}", angle=number * 360 / count)
            for number in range(count)
        ]

    def test_every_label_is_placed(self) -> None:
        placed = place_circular(self._requests(40), CENTRE, 150.0, 11.0)
        assert len(placed.labels) == 40

    def test_no_two_labels_overlap(self) -> None:
        placed = place_circular(self._requests(60), CENTRE, 120.0, 11.0)
        assert placed.collisions() == []

    def test_labels_crowded_at_one_angle_still_do_not_overlap(self) -> None:
        crowded = [Request(text=f"f{number}", angle=90.0) for number in range(25)]
        placed = place_circular(crowded, CENTRE, 150.0, 11.0)
        assert placed.collisions() == []
        assert len(placed.labels) == 25

    def test_a_label_on_the_right_reads_outwards(self) -> None:
        placed = place_circular([Request("x", 90.0)], CENTRE, 100.0, 11.0)
        assert placed.labels[0].side == "right"
        assert placed.labels[0].x > 100.0

    def test_a_label_on_the_left_reads_outwards(self) -> None:
        placed = place_circular([Request("x", 270.0)], CENTRE, 100.0, 11.0)
        assert placed.labels[0].side == "left"
        assert placed.labels[0].box.right <= -100.0

    def test_a_label_with_nothing_in_its_way_is_not_moved(self) -> None:
        placed = place_circular([Request("x", 45.0)], CENTRE, 100.0, 11.0)
        assert placed.labels[0].moved == pytest.approx(0.0)
        assert placed.moved == ()

    def test_how_far_the_furthest_moved_is_reported(self) -> None:
        crowded = [Request(text=f"f{number}", angle=90.0) for number in range(11)]
        placed = place_circular(crowded, CENTRE, 150.0, 11.0)
        assert placed.furthest == pytest.approx(5 * 11.0 * LINE_HEIGHT)

    def test_the_full_text_is_kept_even_when_the_label_is_cut(self) -> None:
        long_name = "gamma-delta insertion target sequence"
        placed = place_circular([Request(long_name, 10.0)], CENTRE, 100.0, 11.0, width=12)
        assert placed.labels[0].text != long_name
        assert placed.labels[0].full == long_name

    def test_the_bounds_hold_every_label(self) -> None:
        placed = place_circular(self._requests(30), CENTRE, 100.0, 11.0)
        bounds = placed.bounds()
        assert bounds is not None
        for label in placed.labels:
            assert bounds.x <= label.box.x
            assert label.box.right <= bounds.right

    def test_no_labels_have_no_bounds(self) -> None:
        assert place_circular([], CENTRE, 100.0, 11.0).bounds() is None

    def test_the_ink_is_shorter_than_the_pitch(self) -> None:
        """Which is why two labels one pitch apart do not touch."""
        assert INK_HEIGHT < LINE_HEIGHT
        placed = place_circular(
            [Request("a", 90.0), Request("b", 91.0)], CENTRE, 150.0, 11.0
        )
        first, second = sorted(placed.labels, key=lambda item: item.y)
        assert first.box.bottom < second.box.y


class TestLinear:
    def _requests(self, count: int, width: float = 800.0) -> list[Along]:
        return [
            Along(text=f"feature {number}", x=number * width / count)
            for number in range(count)
        ]

    def test_every_label_is_placed(self) -> None:
        rows = stack_linear(self._requests(30), top=0.0, font_size=11.0)
        assert len(rows.labels) == 30

    def test_no_two_labels_overlap(self) -> None:
        rows = stack_linear(self._requests(40), top=0.0, font_size=11.0)
        assert rows.collisions() == []

    def test_crowded_labels_are_stacked_into_rows(self) -> None:
        crowded = [Along(text=f"feature {number}", x=400.0) for number in range(8)]
        rows = stack_linear(crowded, top=0.0, font_size=11.0)
        assert rows.rows == 8
        assert rows.collisions() == []

    def test_a_label_stays_over_what_it_names(self) -> None:
        rows = stack_linear([Along(text="EGFP", x=400.0)], top=0.0, font_size=11.0)
        assert rows.labels[0].x == pytest.approx(400.0)
        assert rows.labels[0].anchor == pytest.approx(400.0)

    def test_a_label_at_the_edge_is_pulled_back_inside(self) -> None:
        """A label drawn off the canvas is a label lost."""
        rows = stack_linear(
            [Along(text="right ITR", x=800.0)],
            top=0.0,
            font_size=11.0,
            left=0.0,
            right=800.0,
        )
        assert rows.labels[0].box.right <= 800.0 + 1e-9

    def test_the_leader_target_travels_with_the_label(self) -> None:
        rows = stack_linear(
            [Along(text="x", x=10.0, target=42.0)], top=0.0, font_size=11.0
        )
        assert rows.labels[0].target == 42.0


class TestWidth:
    def test_a_wider_font_estimates_a_wider_label(self) -> None:
        assert estimate_width("abcd", 20.0) > estimate_width("abcd", 10.0)

    def test_the_estimate_is_used_for_placing_and_for_checking(self) -> None:
        """So the two always agree, whatever a renderer actually does."""
        placed = place_circular([Request("abcdefgh", 90.0)], CENTRE, 100.0, 11.0)
        assert placed.labels[0].box.width == pytest.approx(
            estimate_width("abcdefgh", 11.0)
        )


class TestBox:
    def test_boxes_that_touch_do_not_overlap(self) -> None:
        assert not Box(0, 0, 10, 10).overlaps(Box(10, 0, 10, 10))

    def test_padding_makes_touching_boxes_collide(self) -> None:
        assert Box(0, 0, 10, 10).overlaps(Box(10, 0, 10, 10), padding=1.0)
