"""Colours: which kind gets which, and what is refused."""

from __future__ import annotations

import re

import pytest

from plasmid_map.palette import (
    COLOURS,
    GROUPS,
    PALETTES,
    ColourError,
    Palette,
    check_colour,
    group_of,
    is_default,
)


class TestGroups:
    def test_coding_kinds_share_a_group(self) -> None:
        assert group_of("CDS") == group_of("gene") == "coding"

    def test_an_unknown_kind_falls_into_the_default_group(self) -> None:
        assert group_of("something_new") == "other"
        assert is_default("something_new")

    def test_a_known_kind_is_not_a_default(self) -> None:
        assert not is_default("CDS")

    def test_every_group_named_has_a_colour_in_every_palette(self) -> None:
        for name in PALETTES:
            for group in set(GROUPS.values()):
                assert group in COLOURS[name]

    def test_every_palette_has_a_default(self) -> None:
        for name in PALETTES:
            assert "other" in COLOURS[name]


class TestPalette:
    def test_the_two_palettes_differ(self) -> None:
        assert Palette("colour").colour("CDS") != Palette("mono").colour("CDS")

    def test_the_mono_palette_is_grey(self) -> None:
        """Every channel the same, which is what makes it print as grey."""
        for _kind, _group, colour in Palette("mono").rows():
            red, green, blue = (colour[1:3], colour[3:5], colour[5:7])
            assert red == green == blue

    def test_an_override_wins(self) -> None:
        assert Palette("colour", {"CDS": "#123456"}).colour("CDS") == "#123456"

    def test_an_override_leaves_other_kinds_alone(self) -> None:
        palette = Palette("colour", {"CDS": "#123456"})
        assert palette.colour("gene") == COLOURS["colour"]["coding"]

    def test_a_palette_that_does_not_exist_is_refused_by_name(self) -> None:
        with pytest.raises(ColourError) as raised:
            Palette("rainbow")
        assert "colour" in str(raised.value)

    def test_every_kind_is_listed_once(self) -> None:
        rows = Palette().rows()
        assert len(rows) == len(GROUPS)
        assert len({kind for kind, _group, _colour in rows}) == len(GROUPS)

    def test_a_group_colour_can_be_asked_for_directly(self) -> None:
        assert Palette().group_colour("coding") == COLOURS["colour"]["coding"]

    def test_a_group_that_does_not_exist_falls_back(self) -> None:
        assert Palette().group_colour("nonsense") == COLOURS["colour"]["other"]


class TestColours:
    @pytest.mark.parametrize(
        "written", ["#abc", "#AABBCC", "#aabbccdd", "#abcd", "red", "navy"]
    )
    def test_a_colour_that_every_renderer_agrees_on_is_kept(
        self, written: str
    ) -> None:
        assert check_colour(written)

    @pytest.mark.parametrize(
        "written",
        ["blurple", "#12345", "rgb(1,2,3)", "", "#", "url(#x)", "#gggggg"],
    )
    def test_anything_else_is_refused_with_an_example(self, written: str) -> None:
        with pytest.raises(ColourError) as raised:
            check_colour(written)
        assert "#rrggbb" in str(raised.value)

    def test_hexadecimal_is_lowercased_so_a_map_draws_the_same_bytes(self) -> None:
        assert check_colour("#AABBCC") == "#aabbcc"

    def test_every_colour_in_the_table_is_a_colour(self) -> None:
        for name in PALETTES:
            for colour in COLOURS[name].values():
                assert re.fullmatch(r"#[0-9a-f]{6}", colour), colour
