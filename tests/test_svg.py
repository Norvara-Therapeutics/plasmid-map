"""The SVG file: what is in it, and what is deliberately not."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ElementTree

from conftest import construct, spread

from plasmid_map.circular import draw
from plasmid_map.drawing import Style, select
from plasmid_map.svg import (
    FONT,
    document,
    element,
    escape,
    group,
    path,
    rect,
    stylesheet,
    text,
)


def _map() -> str:
    held = construct(features=spread(), name="pMAP1")
    return draw(held, select(held), Style(legend=True)).svg


class TestEscaping:
    def test_angle_brackets_are_escaped(self) -> None:
        assert escape("<script>") == "&lt;script&gt;"

    def test_quotes_are_escaped_so_an_attribute_cannot_be_broken_out_of(
        self,
    ) -> None:
        assert '"' not in escape('a "quoted" name')

    def test_a_label_holding_markup_stays_a_label(self) -> None:
        drawn = text("<b>bold</b>", 0, 0)
        assert "<b>" not in drawn
        assert "&lt;b&gt;" in drawn


class TestElements:
    def test_an_element_with_no_body_closes_itself(self) -> None:
        assert element("rect", {"x": 1}).endswith("/>")

    def test_an_empty_attribute_is_left_out(self) -> None:
        assert "fill" not in element("path", {"d": "M0 0", "fill": ""})

    def test_an_attribute_that_is_none_is_left_out(self) -> None:
        assert "fill" not in element("path", {"d": "M0 0", "fill": None})

    def test_a_path_can_carry_a_tooltip(self) -> None:
        assert "<title>bla 1..10</title>" in path("M0 0", title="bla 1..10")

    def test_a_path_without_a_tooltip_closes_itself(self) -> None:
        assert path("M0 0").endswith("/>")

    def test_a_rectangle_writes_its_numbers_the_one_way(self) -> None:
        assert 'width="10.00"' in rect(0, 0, 10, 10)

    def test_an_empty_group_draws_nothing(self) -> None:
        assert group([]) == ""

    def test_a_group_holds_what_it_is_given(self) -> None:
        assert group(["<rect/>"], "features").startswith('<g class="features"')

    def test_text_carries_its_full_form_as_a_tooltip(self) -> None:
        drawn = text("cut\u2026", 0, 0, title="cut short label")
        assert "<title>cut short label</title>" in drawn

    def test_text_the_same_as_its_title_carries_no_tooltip(self) -> None:
        assert "<title>" not in text("bla", 0, 0, title="bla")


class TestDocument:
    def test_the_document_is_well_formed_xml(self) -> None:
        ElementTree.fromstring(_map().split("\n", 1)[1])

    def test_it_declares_the_svg_namespace(self) -> None:
        assert 'xmlns="http://www.w3.org/2000/svg"' in _map()

    def test_it_carries_a_viewbox_and_a_size(self) -> None:
        drawn = _map()
        assert "viewBox=" in drawn
        assert re.search(r'width="\d', drawn)

    def test_it_carries_a_title_for_a_screen_reader(self) -> None:
        assert "<title>pMAP1, 1,000 bp, circular</title>" in _map()

    def test_the_stylesheet_is_inside_the_file(self) -> None:
        assert "<style>" in _map()

    def test_there_is_no_script(self) -> None:
        assert "<script" not in _map()
        assert "javascript" not in _map().lower()

    def test_nothing_is_fetched_over_the_network(self) -> None:
        drawn = _map()
        assert "http://" not in drawn.replace("http://www.w3.org", "")
        assert "https://" not in drawn
        assert "@import" not in drawn

    def test_the_font_is_named_rather_than_embedded(self) -> None:
        drawn = _map()
        assert FONT in drawn
        assert "@font-face" not in drawn

    def test_the_sizes_asked_for_reach_the_stylesheet(self) -> None:
        sheet = stylesheet(13.0, 8.0, 20.0)
        assert "13.00px" in sheet and "8.00px" in sheet and "20.00px" in sheet

    def test_a_document_with_a_description_carries_it(self) -> None:
        drawn = document([], 10, 10, "t", description="a description")
        assert "<desc>a description</desc>" in drawn

    def test_every_file_ends_with_a_newline(self) -> None:
        assert _map().endswith("</svg>\n")
