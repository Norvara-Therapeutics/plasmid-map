"""The linear feature track.

A horizontal axis with the construct laid along it, the features stacked into
rows above it so that none covers another, their labels stacked above those,
and the ruler and the GC track below. It is the shape to use for a cassette,
an insert, a viral genome or a region of one: anything whose ends are ends
rather than a join.

A circular construct can be drawn this way too, which cuts it at the origin. A
feature crossing that cut is drawn as the two pieces it becomes, at each end
of the track rather than as one block running off the right-hand edge, and the
summary says which features that happened to. The alternative - drawing the
part that runs past the end past the end - would put a feature outside the
construct it is on.
"""

from __future__ import annotations

from typing import Sequence

from plasmid_map.composition import Profile
from plasmid_map.drawing import Drawing, DrawingError, Selection, Style, groups_of
from plasmid_map.feature import Construct, Span, location_text, minimum_width
from plasmid_map.geometry import Box, Point, number, polyline
from plasmid_map.labels import Along, INK_HEIGHT, estimate_width, stack_linear
from plasmid_map.lanes import assign, lane_rows
from plasmid_map.ticks import minor_coordinates, ruler, tick_step
from plasmid_map.svg import document, group, path, rect, text

#: How tall the axis the construct is drawn on is.
AXIS: float = 2.0

#: How much room the ruler takes below the axis.
RULER_HEIGHT: float = 22.0

#: How long a tick is.
TICK: float = 5.0

#: The gap between the labels and the features, and between blocks and rows.
GAP: float = 8.0

#: How much room one row of a legend takes.
LEGEND_ROW: float = 15.0

#: How long an arrowhead is, as a share of the block's thickness.
HEAD_SHARE: float = 0.7

#: How many ticks a linear ruler aims for, before its labels are checked.
TICK_COUNTS: tuple[int, ...] = (10, 8, 6, 5, 4, 3, 2)


def scale(position: float, length: int, left: float, width: float) -> float:
    """Where a base sits along the track, in user units."""
    if length <= 0:
        raise DrawingError("a construct with no length has no scale")
    return left + width * position / length


def pieces_for(
    spans: Sequence[Span], length: int, circular: bool = True
) -> list[tuple[Span, bool]]:
    """The blocks a feature is drawn as, in the order it reads along the DNA.

    Each block comes with whether a join line should be drawn from it to the
    next. Order is the feature's own, not left-to-right on the page, which is
    what puts the arrowhead on the end the feature really ends at: a
    plus-strand gene written ``join(951..1000,1..50)`` ends at base 50, at the
    left-hand edge of the track, and that is where its point belongs.

    A span past the end of a circular construct crosses the origin, and a
    linear track cuts there, so it becomes two blocks with no line between
    them: they are continuous through the cut, not through the middle of the
    sequence, and a dashed line across the middle would say the opposite.

    On a linear construct there is no origin to cross, so a span past the end
    can only be one that was widened to be visible - a one base feature at the
    last base, drawn ten bases wide - and it is clamped to the end instead.
    """
    pieces: list[tuple[Span, bool]] = []
    for span in spans:
        if span.end <= length:
            pieces.append((span, True))
            continue
        if not circular:
            pieces.append((Span(min(span.start, max(length - 1, 0)), length), True))
            continue
        pieces.append((Span(span.start, length), False))
        if span.end - length > 0:
            pieces.append((Span(0, span.end - length), True))
    return pieces


def block(
    left: float,
    right: float,
    top: float,
    thickness: float,
    strand: int,
    colour: str,
    title: str,
) -> str:
    """One feature as a block, pointed at the end it reads towards."""
    head = min(thickness * HEAD_SHARE, max(right - left, 0.0))
    middle = top + thickness / 2
    if strand > 0 and head > 0:
        corner = right - head
        shape = (
            f"M{number(left)} {number(top)}L{number(corner)} {number(top)}"
            f"L{number(right)} {number(middle)}"
            f"L{number(corner)} {number(top + thickness)}"
            f"L{number(left)} {number(top + thickness)}Z"
        )
    elif strand < 0 and head > 0:
        corner = left + head
        shape = (
            f"M{number(right)} {number(top)}L{number(corner)} {number(top)}"
            f"L{number(left)} {number(middle)}"
            f"L{number(corner)} {number(top + thickness)}"
            f"L{number(right)} {number(top + thickness)}Z"
        )
    else:
        shape = (
            f"M{number(left)} {number(top)}L{number(right)} {number(top)}"
            f"L{number(right)} {number(top + thickness)}"
            f"L{number(left)} {number(top + thickness)}Z"
        )
    return path(shape, title=title, fill=colour, **{"class": "feature"})


def _ticks_fit(marks, length: int, left: float, width: float, style: Style) -> bool:
    """Whether the tick labels along the axis keep out of each other's way."""
    from plasmid_map.labels import estimate_width

    boxes: list[Box] = []
    height = style.small_size * INK_HEIGHT
    for tick in marks:
        x = scale(tick.position, length, left, width)
        size = estimate_width(tick.label, style.small_size)
        box = Box(x - size / 2, 0.0, size, height)
        if any(box.overlaps(other, 4.0) for other in boxes):
            return False
        boxes.append(box)
    return True


def choose_ticks(construct: Construct, left: float, width: float, style: Style):
    """The finest ruler whose labels fit along the axis, and its step."""
    for wanted in TICK_COUNTS:
        step = tick_step(construct.length, wanted)
        marks = ruler(construct.length, step, wanted, circular=False)
        if _ticks_fit(marks, construct.length, left, width, style):
            return marks, step
    step = tick_step(construct.length, TICK_COUNTS[-1])
    return ruler(construct.length, step, TICK_COUNTS[-1], circular=False), step


def _gc_track(
    profile: Profile,
    length: int,
    left: float,
    width: float,
    top: float,
    height: float,
) -> list[str]:
    """The GC track as a filled curve about the mean, below the ruler."""
    if profile.highest is None or profile.lowest is None:
        return []
    reach = max((profile.highest - profile.lowest) / 2, 0.01)
    mean = profile.mean or 0.0
    middle = top + height / 2
    points = []
    for window in profile.windows:
        share = window.gc if window.gc is not None else mean
        offset = (share - mean) / reach * (height / 2)
        points.append(
            Point(
                scale(min(window.middle, length), length, left, width),
                middle - max(min(offset, height / 2), -height / 2),
            )
        )
    if not points:
        return []
    area = (
        polyline(points)
        + f"L{number(points[-1].x)} {number(middle)}"
        + f"L{number(points[0].x)} {number(middle)}Z"
    )
    return [
        path(area, **{"class": "gc-fill"}),
        path(
            f"M{number(left)} {number(middle)}L{number(left + width)} {number(middle)}",
            **{"class": "gc-axis"},
        ),
        path(polyline(points), **{"class": "gc-line"}),
    ]


def _skew_track(
    profile: Profile,
    length: int,
    left: float,
    width: float,
    top: float,
    height: float,
) -> list[str]:
    """The GC skew track: one bar per window, coloured by which way it leans.

    Bars rather than a curve, for the same reason as on a circular map: the
    useful thing about skew is its sign and where that changes, and two
    colours say so at a glance.
    """
    limit = profile.skew_limit or 1.0
    middle = top + height / 2
    half = profile.step / 2
    parts = [
        path(
            f"M{number(left)} {number(middle)}L{number(left + width)} "
            f"{number(middle)}",
            **{"class": "gc-axis"},
        )
    ]
    for window in profile.windows:
        if window.skew is None:
            continue
        reach = window.skew / limit * (height / 2)
        if abs(reach) < 0.05:
            continue
        first = scale(max(window.middle - half, 0.0), length, left, width)
        last = scale(min(window.middle + half, float(length)), length, left, width)
        upper, lower = sorted((middle, middle - reach))
        parts.append(
            rect(
                first,
                upper,
                max(last - first, 0.4),
                lower - upper,
                **{"class": "skew-up" if window.skew > 0 else "skew-down"},
            )
        )
    return parts


def draw(
    construct: Construct,
    selection: Selection,
    style: Style,
    profile: Profile | None = None,
) -> Drawing:
    """A linear track of a construct: a standalone SVG and an account of it."""
    if style.start:
        raise DrawingError(
            "--start rotates a circular map and has nothing to turn on a "
            "linear one; draw it with the circular command, or cut the "
            "sequence before drawing it"
        )
    length = construct.length
    minimum = minimum_width(length)
    layout = assign(selection.kept, length, construct.is_circular, minimum)
    lanes = max(layout.lanes, 1)
    left = style.margin
    width = style.width - 2 * style.margin
    if width <= 40:
        raise DrawingError(f"a track {number(style.width)} wide has no room to draw in")

    split = tuple(
        item.feature
        for item in layout.placements
        if item.feature.crosses_origin(length)
    )

    requests: list[Along] = []
    if style.labels != "none":
        for item in layout.placements:
            extent = item.extent
            middle = min((extent.start + extent.end) / 2, length)
            requests.append(
                Along(
                    text=item.feature.label,
                    x=scale(middle, length, left, width),
                    colour=style.palette.colour(item.feature.kind),
                    key=item.feature.kind,
                    target=item.lane * (style.track + style.separation),
                )
            )
    placed = stack_linear(
        requests,
        top=style.margin,
        font_size=style.label_size,
        width=style.label_width,
        left=left,
        right=left + width,
    )
    label_height = placed.rows * placed.row_height if requests else 0.0

    rows_top = style.margin + label_height + (GAP if label_height else 0.0)
    rows = lane_rows(rows_top, style.track, style.separation, lanes)
    axis_y = rows[-1][1] + GAP

    body: list[str] = []
    if requests:
        parts = []
        for label in placed.labels:
            baseline = label.box.bottom - style.label_size * 0.2
            parts.append(
                text(
                    label.text,
                    label.x,
                    baseline,
                    anchor="middle",
                    css="label",
                    title=label.full,
                    fill=label.colour or None,
                )
            )
            parts.append(
                path(
                    f"M{number(label.x)} {number(label.box.bottom + 1)}"
                    f"L{number(label.anchor)} {number(rows_top + label.target - 1)}",
                    **{"class": "leader"},
                )
            )
        body.append(group(parts, "labels"))

    drawn: list[str] = []
    for item in layout.placements:
        top, bottom = rows[item.lane]
        thickness = bottom - top
        pieces = pieces_for(item.spans, length, construct.is_circular)
        where = (
            f"{item.feature.label} ({item.feature.kind}) "
            f"{location_text(item.feature, length)}"
        )
        for index, (span, _bridge) in enumerate(pieces):
            strand = 0
            if item.feature.strand > 0 and index == len(pieces) - 1:
                strand = 1
            elif item.feature.strand < 0 and index == 0:
                strand = -1
            drawn.append(
                block(
                    scale(span.start, length, left, width),
                    scale(span.end, length, left, width),
                    top,
                    thickness,
                    strand,
                    style.palette.colour(item.feature.kind),
                    where,
                )
            )
        for index in range(len(pieces) - 1):
            span, bridge = pieces[index]
            if not bridge:
                continue
            gap_start = span.end
            gap_end = pieces[index + 1][0].start
            if gap_end <= gap_start:
                continue
            middle = top + thickness / 2
            drawn.append(
                path(
                    f"M{number(scale(gap_start, length, left, width))} "
                    f"{number(middle)}"
                    f"L{number(scale(gap_end, length, left, width))} {number(middle)}",
                    **{"class": "join"},
                )
            )
    body.append(group(drawn, "features"))

    body.append(
        group(
            [
                path(
                    f"M{number(left)} {number(axis_y)}"
                    f"L{number(left + width)} {number(axis_y)}",
                    **{"class": "axis"},
                )
            ],
            "axis",
        )
    )

    marks, step = choose_ticks(construct, left, width, style)
    bottom = axis_y
    if style.ruler:
        parts = []
        for coordinate in minor_coordinates(length, step):
            x = scale(coordinate - 1, length, left, width)
            parts.append(
                path(
                    f"M{number(x)} {number(axis_y)}L{number(x)} "
                    f"{number(axis_y + TICK / 2)}",
                    **{"class": "tick-minor"},
                )
            )
        for mark in marks:
            x = scale(mark.position, length, left, width)
            parts.append(
                path(
                    f"M{number(x)} {number(axis_y)}L{number(x)} "
                    f"{number(axis_y + TICK)}",
                    **{"class": "tick"},
                )
            )
            parts.append(
                text(
                    mark.label,
                    x,
                    axis_y + TICK + style.small_size + 1,
                    anchor="middle",
                    css="tick-label",
                )
            )
        body.append(group(parts, "ruler"))
        bottom = axis_y + RULER_HEIGHT

    if (style.gc or style.skew) and profile is None:
        raise DrawingError("a GC track needs the sequence, and this input carries none")
    if style.gc and profile is not None:
        body.append(
            group(
                _gc_track(profile, length, left, width, bottom + GAP, style.gc_height),
                "gc",
            )
        )
        bottom += GAP + style.gc_height
    if style.skew and profile is not None:
        body.append(
            group(
                _skew_track(
                    profile, length, left, width, bottom + GAP, style.gc_height
                ),
                "skew",
            )
        )
        bottom += GAP + style.gc_height

    name = style.title or construct.name
    heading_y = style.margin - 6
    body.append(
        group(
            [
                text(name, left, heading_y, css="title"),
                text(
                    f"{length:,} bp, {construct.topology}",
                    left + estimate_width(name, style.title_size) + 8,
                    heading_y,
                    css="subtitle",
                ),
            ],
            "heading",
        )
    )

    legend = groups_of(layout, style.palette) if style.legend else []
    if legend:
        parts = []
        for position, (group_name, colour) in enumerate(legend):
            y = bottom + GAP + position * LEGEND_ROW
            parts.append(rect(left, y, 10, 10, fill=colour, **{"class": "legend-swatch"}))
            parts.append(text(group_name, left + 14, y + 9, css="legend-label"))
        body.append(group(parts, "legend"))
        bottom += GAP + len(legend) * LEGEND_ROW

    height = bottom + style.margin
    shift = style.title_size + GAP
    svg = document(
        [f'<g transform="translate(0 {number(shift)})">'] + body + ["</g>"],
        width=style.width,
        height=height + shift,
        title=f"{name}, {length:,} bp, {construct.topology}",
        description=construct.description,
        label_size=style.label_size,
        small_size=style.small_size,
        title_size=style.title_size,
    )
    return Drawing(
        svg=svg,
        width=style.width,
        height=height + shift,
        shape="linear",
        construct=construct,
        style=style,
        selection=selection,
        layout=layout,
        placed=placed if requests else None,
        profile=profile if (style.gc or style.skew) else None,
        tick_step=step,
        track=style.track,
        split=split,
    )



