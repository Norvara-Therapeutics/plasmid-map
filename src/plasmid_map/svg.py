"""The SVG file itself: its wrapper, its stylesheet and its escaping.

Everything a drawing needs is inside it. There is no script, no web font and
nothing fetched over the network, so a map opens the same from a disk, a mail
attachment, a journal's submission system and an archive years from now, and
it can be edited afterwards in Illustrator or Inkscape because every feature
is a shape with a class rather than a picture of one.

The font is named rather than embedded. Embedding one would make the file
several hundred kilobytes and would mean redistributing a typeface; naming the
usual sans-serif families means the text is drawn in whatever the reader has,
which is what every figure editor expects to be handed anyway.

Nothing here depends on the clock, on a random number or on the order of a
set, so the same construct drawn twice produces the same bytes and two maps
can be compared with ``diff``.
"""

from __future__ import annotations

from html import escape as _escape
from typing import Iterable, Sequence

from plasmid_map.geometry import number

#: The font families asked for, most preferred first. All are system fonts.
FONT: str = "'Helvetica Neue', Helvetica, Arial, 'Liberation Sans', sans-serif"

#: The stylesheet every drawing carries, so that one is a document rather
#: than a document and a file it needs beside it.
STYLESHEET: str = """
.plasmid-map { font-family: %(font)s; }
.plasmid-map text { fill: #22252a; }
.backbone { fill: none; stroke: #4a4f57; stroke-width: 1.4; }
.backbone-fill { fill: #ffffff; stroke: none; }
.feature { stroke: #ffffff; stroke-width: 0.6; fill-rule: evenodd; }
.feature-edge { fill: none; stroke: #22252a; stroke-width: 0.6; opacity: 0.35; }
.join { fill: none; stroke: #6b7280; stroke-width: 0.9; stroke-dasharray: 3 2; }
.leader { fill: none; stroke: #9aa0a6; stroke-width: 0.7; }
.label { font-size: %(label)spx; }
.label-kind { font-size: %(small)spx; fill: #6b7280; }
.tick { fill: none; stroke: #8b9099; stroke-width: 0.8; }
.tick-minor { fill: none; stroke: #c3c8d0; stroke-width: 0.6; }
.tick-label { font-size: %(small)spx; fill: #6b7280; }
.title { font-size: %(title)spx; font-weight: 600; }
.subtitle { font-size: %(small)spx; fill: #6b7280; }
.gc-fill { fill: #8b9099; opacity: 0.5; stroke: none; fill-rule: evenodd; }
.gc-line { fill: none; stroke: #555b63; stroke-width: 0.7; }
.gc-axis { fill: none; stroke: #c3c8d0; stroke-width: 0.6; stroke-dasharray: 2 3; }
.skew-up { fill: #6f9c5a; opacity: 0.75; stroke: none; }
.skew-down { fill: #a4625f; opacity: 0.75; stroke: none; }
.axis { fill: none; stroke: #4a4f57; stroke-width: 1; }
.legend-swatch { stroke: #ffffff; stroke-width: 0.6; }
.legend-label { font-size: %(small)spx; }
.scale-bar { fill: none; stroke: #22252a; stroke-width: 1.2; }
"""


def escape(text: str) -> str:
    """Text made safe to put inside an SVG element or attribute."""
    return _escape(str(text), quote=True)


def stylesheet(label: float, small: float, title: float) -> str:
    """The stylesheet with the font sizes a drawing was asked for."""
    return STYLESHEET % {
        "font": FONT,
        "label": number(label),
        "small": number(small),
        "title": number(title),
    }


def element(name: str, attributes: dict[str, object], body: str = "") -> str:
    """One element, its attributes written in the order they were given."""
    written = " ".join(
        f'{key}="{escape(value)}"'
        for key, value in attributes.items()
        if value is not None and value != ""
    )
    opened = f"<{name} {written}" if written else f"<{name}"
    if body:
        return f"{opened}>{body}</{name}>"
    return f"{opened}/>"


def path(definition: str, title: str = "", **attributes: object) -> str:
    """A path element, with an optional tooltip carrying its exact extent.

    The tooltip is what lets a map be read precisely without every feature
    having to carry its coordinates in the drawing: a browser shows it on
    hover, and it costs a few bytes rather than a line of text.
    """
    body = f"<title>{escape(title)}</title>" if title else ""
    return element("path", {"d": definition, **attributes}, body)


def rect(
    x: float, y: float, width: float, height: float, **attributes: object
) -> str:
    """A rectangle, for a legend swatch or a linear feature."""
    return element(
        "rect",
        {
            "x": number(x),
            "y": number(y),
            "width": number(width),
            "height": number(height),
            **attributes,
        },
    )


def text(
    body: str,
    x: float,
    y: float,
    anchor: str = "start",
    css: str = "label",
    title: str = "",
    **attributes: object,
) -> str:
    """A piece of text, with an optional tooltip carrying its full form."""
    inner = escape(body)
    if title and title != body:
        inner = f"<title>{escape(title)}</title>{inner}"
    return element(
        "text",
        {
            "class": css,
            "x": number(x),
            "y": number(y),
            "text-anchor": anchor if anchor != "start" else None,
            **attributes,
        },
        inner,
    )


def group(body: Iterable[str], css: str = "", **attributes: object) -> str:
    """A group of elements, which is how a map stays editable afterwards."""
    joined = "".join(body)
    if not joined:
        return ""
    return element("g", {"class": css or None, **attributes}, joined)


def document(
    body: Sequence[str],
    width: float,
    height: float,
    title: str,
    description: str = "",
    label_size: float = 11.0,
    small_size: float = 9.0,
    title_size: float = 15.0,
) -> str:
    """A whole standalone SVG file.

    The viewBox and the width together are what make the file scale: a reader
    can drop it into a two-column layout at any size and the text stays the
    size it was drawn relative to the map.
    """
    header = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'class="plasmid-map" role="img" '
        f'viewBox="0 0 {number(width)} {number(height)}" '
        f'width="{number(width)}" height="{number(height)}" '
        f'aria-label="{escape(title)}">'
    )
    parts = [
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
        header,
        f"<title>{escape(title)}</title>",
    ]
    if description:
        parts.append(f"<desc>{escape(description)}</desc>")
    parts.append(f"<style>{stylesheet(label_size, small_size, title_size)}</style>")
    parts.extend(body)
    parts.append("</svg>")
    return "\n".join(parts) + "\n"
