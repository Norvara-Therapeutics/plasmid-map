"""The standalone HTML view."""

from __future__ import annotations

import base64
import re

from conftest import construct, sequence, spread

from plasmid_map import rows as tabular
from plasmid_map.circular import draw
from plasmid_map.composition import measure
from plasmid_map.drawing import Style, select
from plasmid_map.report import build, data_uri
from plasmid_map.text import notes


def _drawing(**style):
    held = construct(features=spread(), name="pMAP1", sequence=sequence())
    return draw(held, select(held), Style(**style), measure(held.sequence, True))


def _page(**style) -> str:
    return build(_drawing(**style), "x.gb", "genbank")


class TestStandalone:
    def test_it_is_a_whole_document(self) -> None:
        page = _page()
        assert page.startswith("<!DOCTYPE html>")
        assert page.rstrip().endswith("</html>")

    def test_the_drawing_is_inline_rather_than_linked(self) -> None:
        page = _page()
        assert "<svg" in page
        assert "<img" not in page

    def test_the_xml_declaration_is_dropped_when_the_svg_goes_inline(self) -> None:
        assert "<?xml" not in _page()

    def test_the_stylesheet_is_inside_the_page(self) -> None:
        assert "<style>" in _page()

    def test_there_is_no_script(self) -> None:
        page = _page()
        assert "<script" not in page
        assert "javascript" not in page.lower()

    def test_nothing_is_fetched_over_the_network(self) -> None:
        page = _page().replace("http://www.w3.org", "")
        assert "http://" not in page
        assert "https://" not in page
        assert "@import" not in page


class TestContents:
    def test_the_construct_is_named(self) -> None:
        assert "pMAP1" in _page()

    def test_every_note_is_on_the_page(self) -> None:
        drawing = _drawing()
        page = build(drawing, "x.gb", "genbank")
        for note in notes(drawing):
            assert note.replace("'", "&#x27;") in page

    def test_a_map_with_nothing_to_say_says_so(self) -> None:
        from conftest import feature

        held = construct(
            topology="linear", features=(feature(kind="CDS", start=10, end=500),)
        )
        page = build(draw(held, select(held), Style()), "x.gb", "genbank")
        assert "Nothing needed saying" in page

    def test_the_gc_profile_is_offered_when_one_was_drawn(self) -> None:
        assert "GC profile" in _page(gc=True, width=900)

    def test_there_is_no_gc_section_when_no_track_was_drawn(self) -> None:
        assert "GC profile" not in _page()

    def test_a_feature_name_holding_markup_stays_text(self) -> None:
        from conftest import feature

        held = construct(features=(feature(name="<b>bold</b>", start=0, end=100),))
        page = build(draw(held, select(held), Style()), "x.gb", "genbank")
        assert "<b>bold</b>" not in page
        assert "&lt;b&gt;" in page


class TestDownloads:
    def test_a_download_carries_the_table_the_file_would_have_held(self) -> None:
        drawing = _drawing()
        page = build(drawing, "x.gb", "genbank")
        expected = tabular.render(
            tabular.FEATURE_COLUMNS,
            tabular.feature_rows(
                drawing.layout.placements, drawing.construct.length
            ),
        )
        found = re.search(
            r'download="features.tsv" href="data:text/tab-separated-values;'
            r'base64,([A-Za-z0-9+/=]+)"',
            page,
        )
        assert found is not None
        assert base64.b64decode(found.group(1)).decode("utf-8") == expected

    def test_the_gc_profile_downloads_as_the_same_table(self) -> None:
        drawing = _drawing(gc=True, width=900)
        page = build(drawing, "x.gb", "genbank")
        expected = tabular.render(
            tabular.WINDOW_COLUMNS, tabular.window_rows(drawing.profile)
        )
        found = re.search(
            r'download="gc_profile.tsv" href="data:[^;]+;base64,([A-Za-z0-9+/=]+)"',
            page,
        )
        assert found is not None
        assert base64.b64decode(found.group(1)).decode("utf-8") == expected

    def test_a_filename_is_escaped_into_the_link(self) -> None:
        assert 'download="od&quot;d.tsv"' in data_uri("a\tb\n", 'od"d.tsv')


class TestDeterminism:
    def test_the_same_drawing_gives_the_same_page(self) -> None:
        assert _page() == _page()

    def test_nothing_on_the_page_depends_on_the_clock(self) -> None:
        assert not re.search(r"20\d\d-\d\d-\d\d", _page())
