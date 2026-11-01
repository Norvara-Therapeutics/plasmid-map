"""The tables, built once and rendered twice.

Every table the program can write exists here as columns and rows. The
``--output`` file, the tables offered for download inside the HTML report and
the rows printed to standard output are one function's result, so a table
taken out of a report is the table the file would have held.

Coordinates in every table count from 1 and include both ends, which is what
GenBank writes, what a sequence viewer shows and what a person checking a map
against their own records will have in front of them. Inside the program they
are 0-based and half-open; the conversion happens here, once.

A feature that crosses the origin has an ``end`` smaller than its ``start``,
as it does in a GenBank file, and a ``crosses_origin`` column that says so.
Its ``length`` is the number of bases it holds, which for such a feature is
not the difference between its two coordinates and never was.

Values are written unrounded. A share is written in full rather than to one
decimal place; rounding belongs to what a person reads, which is ``text``.
"""

from __future__ import annotations

from typing import Iterable, Sequence

from plasmid_map.composition import Profile
from plasmid_map.drawing import Drawing
from plasmid_map.feature import STRAND_SYMBOLS, Feature
from plasmid_map.feature import location_text as location
from plasmid_map.lanes import Placement
from plasmid_map.palette import GROUPS, Palette, group_of

#: What separates fields. Tab, because a feature name may hold a comma.
SEPARATOR: str = "\t"


def field(value: object) -> str:
    """One value, written for a tab-separated file.

    Floats keep full precision, None becomes an empty field, and a tab or a
    newline inside a name becomes a space so that one row stays one row.
    """
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return str(int(value)) if value.is_integer() else repr(value)
    text = str(value)
    for character in ("\t", "\r", "\n"):
        text = text.replace(character, " ")
    return text


def render(columns: Sequence[str], rows: Iterable[Sequence[object]]) -> str:
    """A whole table as text, with its header line and a trailing newline."""
    lines = [SEPARATOR.join(columns)]
    lines.extend(SEPARATOR.join(field(value) for value in row) for row in rows)
    return "\n".join(lines) + "\n"


def coordinate(position: int, length: int) -> int:
    """A 0-based position as the 1-based coordinate a reader would write.

    A position past the end of a circular construct comes back round, which is
    how the far end of an origin-crossing feature is written.
    """
    return (position % length) + 1 if length else position + 1


FEATURE_COLUMNS: tuple[str, ...] = (
    "kind",
    "name",
    "start",
    "end",
    "length",
    "strand",
    "parts",
    "crosses_origin",
    "partial_ends",
    "lane",
    "drawn_length",
    "widened",
    "location",
)


def partial_ends(feature: Feature) -> str:
    """Which ends of a feature the record marked as running off it.

    One of '', 'start', 'end' or 'both'. A single boolean would say that an
    end was partial without saying which, and the two are different claims
    about where a feature really stops.
    """
    first, last = feature.partial
    if first and last:
        return "both"
    if first:
        return "start"
    if last:
        return "end"
    return ""


def feature_columns(qualifiers: Sequence[str] = ()) -> tuple[str, ...]:
    """The feature columns, with a column for each qualifier asked for."""
    return FEATURE_COLUMNS + tuple(qualifiers)


def ends(feature: Feature, length: int) -> tuple[int, int]:
    """The two coordinates a feature's row reports, counting from 1.

    A feature with no length is a site between two bases, and the two
    coordinates are the bases it lies between: ``102^103`` is reported as 102
    and 103 with a length of 0, rather than as 103 and 102, which would read
    as a feature ending before it starts - the one thing this table uses to
    mean that a feature crosses the origin.
    """
    if feature.length == 0:
        return coordinate(feature.start - 1, length), coordinate(feature.start, length)
    return coordinate(feature.start, length), coordinate(feature.end - 1, length)


def feature_rows(
    placements: Sequence[Placement],
    length: int,
    qualifiers: Sequence[str] = (),
) -> list[list[object]]:
    """One row per drawn feature, with where it is and where it was drawn.

    ``qualifiers`` names qualifiers of the feature to carry across as further
    columns, which is how an extra column of an input table, or a ``/product``
    of a GenBank record, reaches the output rather than only the drawing.
    """
    built: list[list[object]] = []
    for item in placements:
        feature = item.feature
        drawn = sum(span.length for span in item.spans)
        first, last = ends(feature, length)
        row: list[object] = [
            feature.kind,
            feature.name,
            first,
            last,
            feature.length,
            STRAND_SYMBOLS[feature.strand],
            len(feature.parts),
            feature.crosses_origin(length),
            partial_ends(feature),
            item.lane,
            drawn,
            item.widened,
            location(feature, length),
        ]
        row.extend(feature.qualifier(key) for key in qualifiers)
        built.append(row)
    return built


WINDOW_COLUMNS: tuple[str, ...] = (
    "start",
    "end",
    "gc",
    "skew",
    "counted",
    "skipped",
)


def window_rows(profile: Profile) -> list[list[object]]:
    """One row per window of a GC profile."""
    return [
        [
            coordinate(window.start, profile.length),
            coordinate(window.end - 1, profile.length),
            window.gc,
            window.skew,
            window.counted,
            window.skipped,
        ]
        for window in profile.windows
    ]


KIND_COLUMNS: tuple[str, ...] = ("kind", "group", "colour", "features")


def kind_rows(palette: Palette, counts: dict[str, int] | None = None):
    """One row per feature kind in play, and how many of each were drawn.

    'In play' is every kind this program knows a colour for, every kind that
    was drawn, and every kind given a colour of its own with ``--colour`` -
    including one this program has never heard of, which is the only way a
    reader could check that the override took.
    """
    seen = counts or {}
    names = sorted(set(GROUPS) | set(seen) | set(palette.overrides))
    return [
        [kind, group_of(kind), palette.colour(kind), seen.get(kind, 0)]
        for kind in names
    ]


SUMMARY_COLUMNS: tuple[str, ...] = ("measure", "value")


def summary_rows(drawing: Drawing) -> list[list[object]]:
    """The figures a map is described by, as rows rather than as prose."""
    construct = drawing.construct
    rows: list[list[object]] = [
        ["name", construct.name],
        ["accession", construct.accession],
        ["length", construct.length],
        ["topology", construct.topology],
        ["shape", drawing.shape],
        ["features_in_record", len(construct.features)],
        ["features_drawn", drawing.drawn],
        ["features_not_drawn", drawing.selection.dropped],
        ["lanes", drawing.layout.lanes],
        ["crossing_origin", len(drawing.crossing)],
        ["partial_ends", len(drawing.partial)],
        ["widened", len(drawing.layout.widened)],
        ["features_added", drawing.added],
        ["tick_step", drawing.tick_step],
        ["canvas_width", drawing.width],
        ["canvas_height", drawing.height],
    ]
    if drawing.placed is not None and hasattr(drawing.placed, "moved"):
        rows.append(["labels_moved", len(drawing.placed.moved)])
        rows.append(["furthest_label_moved", drawing.placed.furthest])
    if drawing.placed is not None and hasattr(drawing.placed, "rows"):
        rows.append(["label_rows", drawing.placed.rows])
    if drawing.profile is not None:
        rows.append(["gc_mean", drawing.profile.mean])
        rows.append(["gc_windows", len(drawing.profile.windows)])
        rows.append(["gc_measured_round_the_origin", drawing.profile.circular])
        rows.append(["gc_window", drawing.profile.window])
        rows.append(["gc_step", drawing.profile.step])
    return rows


def counts(drawing: Drawing) -> dict[str, int]:
    """How many of each kind were drawn."""
    found: dict[str, int] = {}
    for item in drawing.layout.placements:
        found[item.feature.kind] = found.get(item.feature.kind, 0) + 1
    return found
