"""Numbers and tables as a person reads them, and the notes beside them.

Everything here rounds, abbreviates or aligns, and nothing here is used to
write a file: the files are written by ``rows``, unrounded. The split matters
because a figure shown as 53.7% and a figure written as 0.5374913 are the same
figure for different readers, and only one of them should be the one a later
program parses.

The notes are the other half of what this module is for. A map is a picture,
and a picture cannot say that two of its features were drawn wider than they
are, that a label was moved eighty pixels from the thing it names, or that
eleven features were left out because of a filter. The notes say all of that,
every time, whether or not anybody asked.
"""

from __future__ import annotations

from typing import Sequence

from plasmid_map.drawing import Drawing
from plasmid_map.labels import Placed, Rows
from plasmid_map.palette import Palette, is_default

#: How wide a label may be in a terminal table before it is cut.
CELL_WIDTH: int = 40


def thousands(value: int) -> str:
    """A whole number with thousands separators."""
    return f"{value:,}"


def plural(count: int, singular: str, suffix: str = "s") -> str:
    """A count and its noun, made plural when it needs to be."""
    return f"{thousands(count)} {singular}{'' if count == 1 else suffix}"


def agrees(count: int, one: str, many: str) -> str:
    """The form of a word that agrees with a count.

    'One feature are shorter' is the kind of thing a reader stops reading at,
    and every line here is a line somebody is meant to act on.
    """
    return one if count == 1 else many


def number(value: float | None, places: int = 2, missing: str = "n/a") -> str:
    """A number to a fixed number of places, or a word when there is none."""
    if value is None:
        return missing
    if float(value).is_integer():
        return thousands(int(value))
    return f"{value:,.{places}f}"


def share(value: float | None, places: int = 1, missing: str = "n/a") -> str:
    """A proportion as a percentage.

    A share that is not zero but rounds to zero is shown as '<0.1%' rather
    than as '0.0%', because the two say opposite things about whether
    anything was found at all.
    """
    if value is None:
        return missing
    percentage = value * 100
    if percentage > 0 and round(percentage, places) == 0:
        return f"<{10 ** -places:.{places}f}%"
    if percentage < 100 and round(percentage, places) == 100:
        return f">{100 - 10 ** -places:.{places}f}%"
    return f"{percentage:.{places}f}%"


def bases(count: int) -> str:
    """A length in bases, abbreviated once it is long enough to need it."""
    if count < 1_000:
        return f"{count} bp"
    if count < 1_000_000:
        return f"{count / 1_000:,.1f} kb"
    return f"{count / 1_000_000:,.2f} Mb"


def truncate(text: str, width: int = CELL_WIDTH) -> str:
    """Text cut to a width, with an ellipsis when it was cut."""
    collapsed = " ".join(str(text).split())
    if len(collapsed) <= width:
        return collapsed
    return collapsed[: width - 1] + "\u2026"


def cell(value: object, places: int = 3) -> str:
    """One value as a table cell for a person to read.

    A value that does not exist is blank rather than the word None, and a
    float is rounded: the unrounded figure belongs in the file that ``rows``
    writes, not in a column being read off a terminal.
    """
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return str(int(value)) if value.is_integer() else f"{value:.{places}f}"
    return truncate(str(value))


def table(
    columns: Sequence[str],
    rows: Sequence[Sequence[object]],
    numeric: Sequence[int] = (),
    indent: str = "",
) -> str:
    """A table with its columns aligned, as text.

    Numeric columns are right-aligned and everything else left-aligned, which
    is what makes a column of counts readable at a glance. Width is taken from
    the contents, so nothing is ever cut off mid-number.
    """
    if not rows:
        return ""
    cells = [[cell(value) for value in row] for row in rows]
    widths = [
        max(len(str(columns[position])), *(len(row[position]) for row in cells))
        for position in range(len(columns))
    ]
    lines = []
    header = [
        str(columns[position]).rjust(widths[position])
        if position in numeric
        else str(columns[position]).ljust(widths[position])
        for position in range(len(columns))
    ]
    lines.append(indent + "  ".join(header).rstrip())
    lines.append(indent + "  ".join("-" * width for width in widths))
    for row in cells:
        line = [
            row[position].rjust(widths[position])
            if position in numeric
            else row[position].ljust(widths[position])
            for position in range(len(columns))
        ]
        lines.append(indent + "  ".join(line).rstrip())
    return "\n".join(lines)


def definitions(pairs: Sequence[tuple[str, str]], indent: str = "") -> str:
    """Labelled values, one to a line, with the labels in a column."""
    if not pairs:
        return ""
    width = max(len(label) for label, _value in pairs)
    return "\n".join(
        f"{indent}{label.ljust(width)}  {value}".rstrip() for label, value in pairs
    )


def heading(drawing: Drawing, source: str, source_format: str) -> str:
    """The first line: what was drawn, from what, and in what shape."""
    construct = drawing.construct
    return (
        f"{source} ({source_format}, {construct.topology}, "
        f"{thousands(construct.length)} bp) drawn as a "
        f"{drawing.shape} map"
    )


def _labels_line(drawing: Drawing) -> str | None:
    placed = drawing.placed
    if placed is None:
        return None
    if isinstance(placed, Placed):
        moved = len(placed.moved)
        if not moved:
            return f"{plural(len(placed.labels), 'label')}, none moved"
        return (
            f"{plural(len(placed.labels), 'label')}, {moved} moved apart, "
            f"furthest {number(placed.furthest, 0)} px"
        )
    if isinstance(placed, Rows):
        return f"{plural(len(placed.labels), 'label')} in {plural(placed.rows, 'row')}"
    return None


def measures(drawing: Drawing) -> list[tuple[str, str]]:
    """The figures that describe a map, labelled for a person."""
    construct = drawing.construct
    layout = drawing.layout
    pairs: list[tuple[str, str]] = [
        ("features in the record", thousands(len(construct.features))),
        (
            "drawn",
            f"{thousands(drawing.drawn)} in {plural(layout.lanes, 'lane')}"
            if drawing.shape == "circular"
            else f"{thousands(drawing.drawn)} in {plural(layout.lanes, 'row')}",
        ),
        ("not drawn", thousands(drawing.selection.dropped)),
    ]
    if construct.is_circular:
        pairs.append(("crossing the origin", thousands(len(drawing.crossing))))
    if drawing.partial:
        pairs.append(("with an end off the record", thousands(len(drawing.partial))))
    widened = len(layout.widened)
    if widened:
        pairs.append(
            (
                "drawn wider than they are",
                f"{thousands(widened)}, the narrowest drawn as "
                f"{bases(layout.minimum)}",
            )
        )
    line = _labels_line(drawing)
    if line:
        pairs.append(("labels", line))
    if drawing.style.ruler:
        pairs.append(("ruler", f"a tick every {bases(drawing.tick_step)}"))
    pairs.append(
        ("canvas", f"{number(drawing.width, 0)} x {number(drawing.height, 0)} px")
    )
    profile = drawing.profile
    if profile is not None:
        if drawing.style.gc:
            pairs.append(
                (
                    "GC content",
                    f"{share(profile.mean)} overall, {share(profile.lowest)} to "
                    f"{share(profile.highest)} across the windows",
                )
            )
        if drawing.style.skew:
            pairs.append(
                ("GC skew", f"reaches {number(profile.skew_limit, 2)} either way")
            )
        pairs.append(("windows", _windows_line(profile)))
    return pairs


def _windows_line(profile) -> str:
    """How the windows were measured, which is not the same on the two shapes.

    A circular sequence is measured round the origin, so every window is the
    full width and the track closes. A linear one stops, and the last window
    is short. A GC track with a seam at the origin is the usual sign that a
    circular sequence was measured as a linear one, so which was done is said
    rather than left to be inferred from the picture.
    """
    shape = "round the origin" if profile.circular else "stopping at the end"
    tail = ""
    if profile.windows:
        last = profile.windows[-1]
        held = last.end - last.start
        if not profile.circular and held < profile.window:
            tail = f", the last {bases(held)}"
    return (
        f"{plural(len(profile.windows), 'window')} of {bases(profile.window)} "
        f"every {bases(profile.step)}, {shape}{tail}"
    )


def kinds_drawn(drawing: Drawing) -> tuple[list[str], list[list[object]]]:
    """A table of what kinds were drawn, most common first."""
    found: dict[str, int] = {}
    for item in drawing.layout.placements:
        found[item.feature.kind] = found.get(item.feature.kind, 0) + 1
    ordered = sorted(found.items(), key=lambda pair: (-pair[1], pair[0]))
    palette = drawing.style.palette
    return (
        ["kind", "drawn", "colour"],
        [[kind, count, palette.colour(kind)] for kind, count in ordered],
    )


def notes(drawing: Drawing, loaded_records: int = 1, stated: bool = True) -> list[str]:
    """Everything a reader should know that the picture cannot say.

    Each is a fact about this drawing rather than a warning about drawings in
    general, and each names the number it is about, so that a note can be
    checked against the map.
    """
    found: list[str] = []
    construct = drawing.construct
    selection = drawing.selection
    if loaded_records > 1:
        found.append(
            f"The file holds {thousands(loaded_records)} records and the first "
            f"of them, {construct.name}, was drawn. The others were not."
        )
    if not stated:
        found.append(
            f"The input does not record whether the ends are joined, so it was "
            f"drawn as {construct.topology}. Say which it is with --topology."
        )
    if selection.whole_record:
        found.append(
            f"{plural(len(selection.whole_record), 'feature')} describing the "
            f"whole record {agrees(len(selection.whole_record), 'was', 'were')} "
            "not drawn, since a band round the entire construct says nothing "
            "and costs a lane. Draw them with --kind source."
        )
    if selection.hidden:
        listed = ", ".join(
            f"{kind} ({thousands(count)})" for kind, count in selection.hidden
        )
        found.append(f"Left out by the kind filter: {listed}.")
    if selection.short:
        found.append(
            f"{plural(len(selection.short), 'feature')} "
            f"{agrees(len(selection.short), 'was', 'were')} shorter than "
            f"--minimum-length and {agrees(len(selection.short), 'was', 'were')} "
            "left out."
        )
    if drawing.crossing:
        names = ", ".join(sorted({item.label for item in drawing.crossing}))
        if drawing.shape == "circular":
            found.append(
                f"{plural(len(drawing.crossing), 'feature')} "
                f"{agrees(len(drawing.crossing), 'crosses', 'cross')} the origin "
                f"and {agrees(len(drawing.crossing), 'is', 'are')} each drawn as "
                f"one arc through it, not two: {names}. "
                f"{agrees(len(drawing.crossing), 'Its length is', 'Their lengths are')} "
                "the sum of their parts."
            )
        else:
            found.append(
                f"{plural(len(drawing.crossing), 'feature')} "
                f"{agrees(len(drawing.crossing), 'crosses', 'cross')} the origin, "
                f"and a linear track cuts there, so each is drawn as the two "
                f"pieces it becomes, at either end: {names}."
            )
    if drawing.partial:
        found.append(
            f"{plural(len(drawing.partial), 'feature')} "
            f"{agrees(len(drawing.partial), 'has', 'have')} an end the record "
            f"marks as running past what was sequenced; "
            f"{agrees(len(drawing.partial), 'it is', 'they are')} drawn stopping "
            "at the coordinate given."
        )
    widened = drawing.layout.widened
    if widened:
        found.append(
            f"{plural(len(widened), 'feature')} "
            f"{agrees(len(widened), 'is', 'are')} shorter than "
            f"{bases(drawing.layout.minimum)} and {agrees(len(widened), 'is', 'are')} "
            f"drawn that wide so {agrees(len(widened), 'it', 'they')} can be seen, "
            f"so {agrees(len(widened), 'its', 'their')} width on the map is not "
            f"{agrees(len(widened), 'its', 'their')} length."
        )
    placed = drawing.placed
    if isinstance(placed, Placed) and placed.moved:
        found.append(
            f"{plural(len(placed.moved), 'label')} "
            f"{agrees(len(placed.moved), 'was', 'were')} moved along "
            f"{agrees(len(placed.moved), 'its', 'their')} side of the map to keep "
            f"{agrees(len(placed.moved), 'it clear of the rest', 'them apart')}, "
            f"the furthest by {number(placed.furthest, 0)} px. Each is joined to "
            "what it names by a leader line; none was dropped."
        )
    if isinstance(placed, Rows) and placed.rows > 1:
        found.append(
            f"The labels needed {plural(placed.rows, 'row')} to keep them "
            "apart; each is joined to what it names by a leader line."
        )
    if drawing.thinned:
        found.append(
            f"{plural(drawing.layout.lanes, 'lane')} did not fit at the full "
            f"band thickness, so the bands were thinned to "
            f"{number(drawing.track, 1)} px to leave the middle of the map "
            "readable. Widen it with --width to get them back."
        )
    overridden = drawing.style.palette.overrides
    unknown = [
        kind
        for kind in selection.unknown_kinds
        if is_default(kind) and kind not in overridden
    ]
    if unknown:
        found.append(
            f"{agrees(len(unknown), 'This feature kind has', 'These feature kinds have')} "
            f"no colour of their own and "
            f"{agrees(len(unknown), 'was', 'were')} drawn in the default: "
            f"{', '.join(unknown)}."
        )
    if drawing.added:
        found.append(
            f"{plural(drawing.added, 'feature')} "
            f"{agrees(drawing.added, 'was', 'were')} added to the record before "
            "drawing, from a feature table or the command line."
        )
    if isinstance(placed, Rows) and placed.shortened:
        found.append(
            f"{plural(placed.shortened, 'label')} "
            f"{agrees(placed.shortened, 'was', 'were')} wider than the track and "
            f"{agrees(placed.shortened, 'was', 'were')} cut further to stay on "
            "the canvas; the whole text is still in the drawing as a tooltip."
        )
    profile = drawing.profile
    if profile is not None and profile.skipped:
        found.append(
            f"{plural(profile.skipped, 'base')} "
            f"{agrees(profile.skipped, 'is', 'are')} not A, C, G, T or U and "
            f"{agrees(profile.skipped, 'was', 'were')} left out of the GC share "
            "rather than counted against it."
        )
    collisions = drawing.label_collisions
    if collisions:
        found.append(
            f"{plural(len(collisions), 'pair')} of labels still "
            f"{agrees(len(collisions), 'overlaps', 'overlap')}, which should not "
            "happen; please report it."
        )
    return found


def summary(
    drawing: Drawing,
    source: str,
    source_format: str,
    loaded_records: int = 1,
    stated: bool = True,
) -> str:
    """The whole of what goes to standard output after a map is drawn."""
    parts = [heading(drawing, source, source_format)]
    parts.append(definitions(measures(drawing), indent="  "))
    columns, rows = kinds_drawn(drawing)
    if rows:
        parts.append("")
        parts.append(table(columns, rows, numeric=(1,), indent="  "))
    for note in notes(drawing, loaded_records, stated):
        parts.append("")
        parts.append(f"  note: {note}")
    return "\n".join(parts) + "\n"


def palette_table(palette: Palette, counts: dict[str, int] | None = None) -> str:
    """Every feature kind this program knows and the colour it is drawn in."""
    from plasmid_map.rows import KIND_COLUMNS, kind_rows

    return table(
        list(KIND_COLUMNS), kind_rows(palette, counts), numeric=(3,), indent="  "
    )
