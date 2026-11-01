"""Tick spacing and tick labels."""

from __future__ import annotations

import pytest

from plasmid_map.ticks import (
    STEPS,
    RulerError,
    coordinates,
    minor_coordinates,
    ruler,
    tick_label,
    tick_step,
)


class TestStep:
    def test_a_step_is_always_on_the_one_two_five_scale(self) -> None:
        for length in range(50, 200_000, 331):
            step = tick_step(length)
            mantissa = step
            while mantissa % 10 == 0 and mantissa >= 10:
                mantissa //= 10
            assert mantissa in STEPS, (length, step)

    def test_the_step_gives_at_least_as_many_ticks_as_were_asked_for(self) -> None:
        for length in range(50, 200_000, 331):
            assert length / tick_step(length, 10) >= 10 or tick_step(length, 10) == 1

    def test_a_plasmid_of_2686_bases_is_ticked_every_200(self) -> None:
        assert tick_step(2686, 10) == 200

    def test_asking_for_fewer_ticks_gives_a_coarser_step(self) -> None:
        assert tick_step(4361, 4) >= tick_step(4361, 10)

    def test_a_very_short_construct_is_ticked_every_base(self) -> None:
        assert tick_step(5, 10) == 1

    def test_a_construct_with_no_length_has_no_ruler(self) -> None:
        with pytest.raises(RulerError):
            tick_step(0)

    def test_a_ruler_needs_at_least_one_tick(self) -> None:
        with pytest.raises(RulerError):
            tick_step(100, 0)


class TestCoordinates:
    def test_base_one_is_always_marked(self) -> None:
        assert coordinates(4361, 500)[0] == 1

    def test_the_marks_are_multiples_of_the_step(self) -> None:
        marked = coordinates(4361, 500)
        assert marked[1:] == [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]

    def test_a_mark_on_the_last_base_of_a_circle_is_left_out(self) -> None:
        """It is the same place as base 1 and would print two labels on one."""
        assert 1000 not in coordinates(1000, 100, circular=True)
        assert 1000 in coordinates(1000, 100, circular=False)

    def test_a_step_that_does_not_move_along_is_refused(self) -> None:
        with pytest.raises(RulerError):
            coordinates(100, 0)

    def test_minor_marks_fall_between_the_major_ones(self) -> None:
        minor = minor_coordinates(1000, 100)
        assert all(position % 100 for position in minor)
        assert minor[0] == 20

    def test_there_are_no_minor_marks_without_divisions(self) -> None:
        assert minor_coordinates(1000, 100, divisions=1) == []


class TestLabels:
    def test_base_one_keeps_its_number(self) -> None:
        assert tick_label(1, 1000) == "1"

    def test_a_round_thousand_becomes_kilobases(self) -> None:
        assert tick_label(2000, 1000) == "2 kb"

    def test_a_fine_step_keeps_the_number(self) -> None:
        assert tick_label(2000, 200) == "2,000"

    def test_a_number_that_is_not_a_round_thousand_keeps_its_digits(self) -> None:
        assert tick_label(2500, 1000) == "2,500"


class TestRuler:
    def test_a_tick_is_drawn_one_base_before_its_coordinate(self) -> None:
        """The program counts from 0 and a reader counts from 1."""
        ticks = ruler(4361, 500)
        assert ticks[0].position == 0 and ticks[0].coordinate == 1
        assert ticks[1].position == 499 and ticks[1].coordinate == 500

    def test_every_tick_carries_a_label(self) -> None:
        assert all(tick.label for tick in ruler(4361))

    def test_the_number_of_ticks_follows_the_step(self) -> None:
        assert len(ruler(1000, 100, circular=True)) == 10
