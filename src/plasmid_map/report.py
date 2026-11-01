"""The standalone HTML view: the map, what it says, and the tables behind it.

One file with nothing outside it. The drawing is the same SVG that ``--svg``
writes, placed inline rather than linked, the stylesheet is embedded, there is
no script and nothing is fetched, so the page opens the same from a disk, an
attachment or an archive. Each table on it downloads as the tab-separated file
``--output`` would have written, carried as a data URI and produced by the
same serialiser, so the two cannot drift apart.

The same input produces the same bytes. Nothing here reads the clock or the
name of the machine, which is what lets two reports be compared with ``diff``.
"""

from __future__ import annotations

import base64
from html import escape
from typing import Sequence

from plasmid_map import rows as tabular
from plasmid_map import text as prose
from plasmid_map.drawing import Drawing

#: The stylesheet the page carries. Print rules included, because the usual
#: reason to open one of these is to put the map in something else.
STYLESHEET: str = """
:root { color-scheme: light; }
body {
  margin: 0; padding: 2rem 1.5rem 4rem;
  font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
  color: #22252a; background: #ffffff; line-height: 1.5;
}
main { max-width: 68rem; margin: 0 auto; }
h1 { font-size: 1.5rem; margin: 0 0 0.25rem; }
h2 { font-size: 1rem; margin: 2rem 0 0.5rem; font-weight: 600; }
p.lede { margin: 0 0 1.5rem; color: #5c636e; }
figure { margin: 0 0 1.5rem; }
figure svg { max-width: 100%; height: auto; display: block; }
dl.measures { display: grid; grid-template-columns: auto 1fr; gap: 0.2rem 1.25rem;
  margin: 0; }
dl.measures dt { color: #5c636e; }
dl.measures dd { margin: 0; }
ul.notes { margin: 0; padding-left: 1.1rem; color: #3b4048; }
ul.notes li { margin-bottom: 0.45rem; }
table { border-collapse: collapse; font-size: 0.85rem; width: 100%; }
th, td { text-align: left; padding: 0.3rem 0.6rem; border-bottom: 1px solid #e4e7ec; }
th { font-weight: 600; color: #5c636e; }
td.number, th.number { text-align: right; font-variant-numeric: tabular-nums; }
.swatch { display: inline-block; width: 0.7rem; height: 0.7rem;
  border: 1px solid #ffffff; vertical-align: -1px; margin-right: 0.35rem; }
a.download { font-size: 0.8rem; color: #3d6fb4; text-decoration: none; }
a.download:hover { text-decoration: underline; }
footer { margin-top: 3rem; color: #8b9099; font-size: 0.8rem; }
@media print {
  body { padding: 0; }
  a.download { display: none; }
}
"""


def data_uri(body: str, name: str) -> str:
    """A table as a link that downloads it, carried inside the page."""
    encoded = base64.b64encode(body.encode("utf-8")).decode("ascii")
    return (
        f'<a class="download" download="{escape(name)}" '
        f'href="data:text/tab-separated-values;base64,{encoded}">'
        f"download {escape(name)}</a>"
    )


def _table(
    columns: Sequence[str],
    body: Sequence[Sequence[object]],
    numeric: Sequence[int] = (),
    swatch: int | None = None,
) -> str:
    """One table, with the same values the tab-separated file would hold."""
    right = ' class="number"'
    head = "".join(
        "<th" + (right if index in numeric else "") + f">{escape(name)}</th>"
        for index, name in enumerate(columns)
    )
    lines = []
    for row in body:
        cells = []
        for index, value in enumerate(row):
            shown = escape(prose.cell(value))
            if index == swatch and isinstance(value, str) and value.startswith("#"):
                shown = (
                    f'<span class="swatch" style="background:{escape(value)}"></span>'
                    f"{shown}"
                )
            cells.append(
                "<td" + (right if index in numeric else "") + f">{shown}</td>"
            )
        lines.append("<tr>" + "".join(cells) + "</tr>")
    return (
        f"<table><thead><tr>{head}</tr></thead>"
        f"<tbody>{''.join(lines)}</tbody></table>"
    )


def build(
    drawing: Drawing,
    source: str,
    source_format: str,
    loaded_records: int = 1,
    stated: bool = True,
) -> str:
    """The whole report as one HTML document."""
    construct = drawing.construct
    name = drawing.style.title or construct.name
    title = f"{name} - plasmid map"
    inline = drawing.svg
    if inline.startswith("<?xml"):
        inline = inline.split("\n", 1)[1]

    measures = "".join(
        f"<dt>{escape(label)}</dt><dd>{escape(value)}</dd>"
        for label, value in prose.measures(drawing)
    )
    told = prose.notes(drawing, loaded_records, stated)
    notes = (
        "<ul class='notes'>"
        + "".join(f"<li>{escape(note)}</li>" for note in told)
        + "</ul>"
        if told
        else "<p class='lede'>Nothing needed saying about this one.</p>"
    )

    feature_body = tabular.feature_rows(drawing.layout.placements, construct.length)
    feature_table = tabular.render(tabular.FEATURE_COLUMNS, feature_body)
    kind_columns, kind_body = prose.kinds_drawn(drawing)

    sections = [
        "<h2>What was drawn</h2>",
        f"<dl class='measures'>{measures}</dl>",
        "<h2>Notes</h2>",
        notes,
        "<h2>Feature kinds</h2>",
        _table(kind_columns, kind_body, numeric=(1,), swatch=2),
        "<h2>Features</h2>",
        f"<p class='lede'>{data_uri(feature_table, 'features.tsv')}</p>",
        _table(
            list(tabular.FEATURE_COLUMNS),
            feature_body,
            numeric=(2, 3, 4, 6, 9, 10),
        ),
    ]
    if drawing.profile is not None:
        window_body = tabular.window_rows(drawing.profile)
        described = (
            f"{len(window_body):,} windows of "
            f"{prose.bases(drawing.profile.window)}, every "
            f"{prose.bases(drawing.profile.step)}."
        )
        link = data_uri(
            tabular.render(tabular.WINDOW_COLUMNS, window_body), "gc_profile.tsv"
        )
        sections.extend(
            [
                "<h2>GC profile</h2>",
                f"<p class='lede'>{escape(described)} {link}</p>",
            ]
        )

    return (
        "<!DOCTYPE html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{escape(title)}</title>\n"
        f"<style>{STYLESHEET}</style>\n</head>\n<body>\n<main>\n"
        f"<h1>{escape(name)}</h1>\n"
        f"<p class='lede'>{escape(construct.description or source)}</p>\n"
        f"<figure>{inline}</figure>\n"
        + "\n".join(sections)
        + "\n<footer>Drawn by plasmid-map from "
        f"{escape(source)} ({escape(source_format)}). "
        "Every figure on this page is in the tables above it.</footer>\n"
        "</main>\n</body>\n</html>\n"
    )
