"""GC content and skew, and the windows they are measured in."""

from __future__ import annotations

import pytest
from conftest import sequence

from plasmid_map.composition import (
    CompositionError,
    at_rich,
    default_step,
    default_window,
    gc_share,
    gc_skew,
    measure,
    shares,
)


class TestShare:
    def test_a_half_gc_sequence_is_a_half(self) -> None:
        share, counted, skipped = gc_share("ACGT")
        assert share == 0.5
        assert counted == 4
        assert skipped == 0

    def test_ambiguity_is_left_out_of_both_sides(self) -> None:
        """Counting N against GC would make a run of N look AT-rich."""
        share, counted, skipped = gc_share("GCNNNN")
        assert share == 1.0
        assert counted == 2
        assert skipped == 4

    def test_a_stretch_with_nothing_countable_has_no_share(self) -> None:
        share, counted, _skipped = gc_share("NNNN")
        assert share is None
        assert counted == 0

    def test_case_does_not_matter(self) -> None:
        assert gc_share("acgt")[0] == gc_share("ACGT")[0]

    def test_uracil_counts_as_a_base(self) -> None:
        assert gc_share("ACGU")[1] == 4


class TestSkew:
    def test_more_g_than_c_leans_positive(self) -> None:
        assert gc_skew("GGGC") == 0.5

    def test_an_even_split_is_zero(self) -> None:
        assert gc_skew("GCGC") == 0.0

    def test_neither_g_nor_c_has_no_skew_rather_than_zero(self) -> None:
        """No bases to divide and an even split are different facts."""
        assert gc_skew("ATATAT") is None


class TestWindows:
    def test_a_circular_sequence_has_no_seam(self) -> None:
        profile = measure("ACGT" * 25, circular=True, window=30, step=30)
        assert profile.windows[-1].end > profile.length
        assert all(
            window.end - window.start == 30 for window in profile.windows
        )

    def test_a_linear_sequence_keeps_its_short_last_window(self) -> None:
        """Every base belongs to a window; the final one is simply shorter."""
        profile = measure("ACGT" * 25, circular=False, window=30, step=30)
        assert profile.windows[-1].end == 100
        assert profile.windows[-1].end - profile.windows[-1].start == 10

    def test_every_base_of_a_linear_sequence_is_in_a_window(self) -> None:
        profile = measure(sequence(997), circular=False)
        covered = set()
        for window in profile.windows:
            covered.update(range(window.start, window.end))
        assert covered == set(range(997))

    def test_the_mean_is_taken_over_the_whole_sequence(self) -> None:
        body = sequence(500)
        profile = measure(body, circular=True)
        assert profile.mean == pytest.approx(gc_share(body)[0])

    def test_a_window_crossing_the_origin_takes_bases_from_both_ends(self) -> None:
        profile = measure("G" * 10 + "A" * 80 + "C" * 10, True, 20, 30)
        last = profile.windows[-1]
        assert (last.start, last.end) == (90, 110)
        assert last.gc == 1.0

    def test_the_spread_is_reported(self) -> None:
        profile = measure(sequence(2000), circular=True)
        assert profile.lowest is not None and profile.highest is not None
        assert profile.lowest <= profile.mean <= profile.highest

    def test_the_at_richest_window_can_be_found(self) -> None:
        profile = measure(sequence(2000), circular=True)
        poorest = at_rich(profile, 1)[0]
        assert poorest.gc == profile.lowest

    def test_shares_leaves_out_the_windows_that_have_none(self) -> None:
        profile = measure("N" * 50 + "GCGC" * 50, circular=False, window=25, step=25)
        assert len(shares(profile.windows)) < len(profile.windows)

    def test_a_window_wider_than_the_sequence_is_cut_to_it(self) -> None:
        profile = measure("ACGTACGTAC", circular=False, window=500)
        assert profile.window == 10
        assert len(profile.windows) == 1


class TestDefaults:
    def test_the_window_is_never_finer_than_fifty_bases(self) -> None:
        assert default_window(400) == 50

    def test_a_long_construct_gets_a_proportional_window(self) -> None:
        assert default_window(40_000) == 1_000

    def test_the_step_is_a_fifth_of_the_window(self) -> None:
        assert default_step(100, 10_000) == 20

    def test_the_step_is_coarsened_so_a_profile_stays_a_usable_size(self) -> None:
        assert 10_000_000 // default_step(50, 10_000_000) <= 2000


class TestRefusals:
    def test_an_empty_sequence_is_refused(self) -> None:
        with pytest.raises(CompositionError):
            measure("", circular=True)

    def test_a_window_of_no_bases_is_refused(self) -> None:
        with pytest.raises(CompositionError):
            measure("ACGT", circular=True, window=0)

    def test_a_step_that_does_not_move_along_is_refused(self) -> None:
        with pytest.raises(CompositionError):
            measure("ACGT", circular=True, window=2, step=0)


class TestRealRecords:
    def test_pbr322_has_the_gc_content_it_is_known_for(self, pbr322) -> None:
        profile = measure(pbr322.sequence, pbr322.is_circular)
        assert profile.mean == pytest.approx(0.538, abs=0.002)

    def test_the_aav2_repeats_are_the_gc_richest_part_of_it(self, aav2) -> None:
        """The inverted terminal repeats are hairpins, and hairpins are GC."""
        profile = measure(aav2.sequence, aav2.is_circular, window=145, step=145)
        richest = max(profile.windows, key=lambda window: window.gc)
        assert richest.start < 145 or richest.end > aav2.length - 145
