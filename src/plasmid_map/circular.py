"""The circular map.

Rings from the outside in: the labels, then the feature lanes with the longest
features outermost, then the backbone circle the lanes sit on, then the ruler
ticking inwards, then the GC and skew tracks if they were asked for, and the
construct's name in the middle.

The radius is worked out from the width rather than the other way round. The
labels are measured first - their angles do not depend on the radius, so which
side each falls on and how wide the widest on each side is are both known
before anything is drawn - and the circle is made as large as what is left
allows. A map therefore comes out at the width that was asked for, instead of
at the width plus however long the longest gene name happened to be.

The ruler is coarsened until its own labels stop colliding, so a short
construct is not given a ring of numbers printed over each other.
"""

from __future__ import annotations

from dataclasses import dataclass

from plasmid_map.composition import Profile
from plasmid_map.drawing import (
    Drawing,
    DrawingError,
    Selection,
    Style,
    groups_of,
)
from plasmid_map.feature import Construct, location_text, minimum_width
from plasmid_map.geometry import (
    Box,
    Point,
    angle_at,
    arc_band,
    arc_line,
    number,
    polar,
    polyline,
    radial,
    ring,
    unwrapped_angles,
)
from plasmid_map.labels import (
    LINE_HEIGHT,
    Request,
    estimate_width,
    place_circular,
    shorten,
)
from plasmid_map.lanes import Layout, assign, lane_radii
from plasmid_map.ticks import minor_coordinates, ruler, tick_step
from plasmid_map.svg import document, group, path, rect, text

#: How far a leader line stands off the outside of the map.
LEADER_GAP: float = 14.0

#: How long a major tick is, and how far its label sits inside it.
TICK: float = 6.0
TICK_LABEL_GAP: float = 5.0

#: The numbers of ticks tried, most first, until the tick labels stop
#: colliding with one another.
TICK_COUNTS: tuple[int, ...] = (10, 8, 6, 5, 4, 3, 2)

#: How much room one row of a legend takes.
LEGEND_ROW: float = 15.0

#: The smallest backbone radius a map is drawn at before it is refused.
LEAST_RADIUS: float = 40.0

#: The smallest radius kept clear in the middle for the construct's name,
#: whatever the name is. A longer name reserves more.
LEAST_CENTRE: float = 28.0


def centre_reserve(name: str, subtitle: str, style: Style) -> float:
    """How much of the middle is kept clear for the name and the length.

    A track drawn through the name would make both unreadable, so the two
    lines written in the middle decide how far in the innermost track may
    come, and a track that does not fit is refused rather than drawn over
    them.
    """
    return max(
        LEAST_CENTRE,
        estimate_width(name, style.title_size) / 2 + 8,
        estimate_width(subtitle, style.small_size) / 2 + 8,
    )

#: The longest an arrowhead may be, in degrees, however thick the band is.
LONGEST_HEAD: float = 28.0

#: Degrees in a radian, for turning a band's thickness into an arrowhead.
PER_RADIAN: float = 57.29577951308232


#: The largest share of the radius the feature lanes may take between them.
#: Annotation that needs eight lanes would otherwise leave nothing in the
#: middle for the name, the ruler and the GC track, so the bands are thinned
#: instead and the summary says they were.
BAND_SHARE: float = 0.55

#: How thin a feature band may be made before thinning stops helping.
LEAST_TRACK: float = 3.0


def fit_bands(radius: float, lanes: int, style: Style) -> tuple[float, float]:
    """The band thickness and spacing that leave the middle of the map usable.

    A map of a well-annotated plasmid may need eight or ten lanes. Drawn at
    the full band thickness those lanes would reach the centre, so when they
    would take more than their share of the radius every band and every gap is
    scaled down together, which keeps the proportions of the drawing and
    leaves room for the ruler and the name.
    """
    if lanes <= 0:
        return style.track, style.separation
    needed = (lanes - 1) * (style.track + style.separation) + style.track
    allowed = radius * BAND_SHARE
    if needed <= allowed or needed <= 0:
        return style.track, style.separation
    scale = allowed / needed
    track = max(style.track * scale, LEAST_TRACK)
    return track, style.separation * scale


@dataclass(frozen=True)
class Mark:
    """One tick of the ruler, with the angle it is drawn at."""

    position: int
    label: str
    angle: float


def label_requests(
    construct: Construct, layout: Layout, style: Style
) -> list[Request]:
    """One label request per drawn feature, aimed at the middle of its extent.

    The middle of the extent rather than of the parts: a label for a spliced
    feature belongs opposite the whole of it, not opposite its first exon.
    """
    requests: list[Request] = []
    for item in layout.placements:
        extent = item.extent
        middle = (extent.start + extent.end) / 2
        requests.append(
            Request(
                text=item.feature.label,
                angle=angle_at(middle, construct.length, style.start),
                colour=style.palette.colour(item.feature.kind),
                key=item.feature.kind,
            )
        )
    return requests


def side_room(requests: list[Request], style: Style) -> tuple[float, float]:
    """How much width the labels need on the left and on the right."""
    left = 0.0
    right = 0.0
    for request in requests:
        width = estimate_width(
            shorten(request.text, style.label_width), style.label_size
        )
        if request.angle % 360 < 180:
            right = max(right, width)
        else:
            left = max(left, width)
    return left, right


def _ticks_fit(
    marks: list[Mark], centre: Point, radius: float, style: Style
) -> bool:
    """Whether tick labels at this radius keep out of each other's way."""
    boxes: list[Box] = []
    height = style.small_size * LINE_HEIGHT
    for mark in marks:
        at = polar(centre, radius, mark.angle)
        width = estimate_width(mark.label, style.small_size)
        box = Box(at.x - width / 2, at.y - height / 2, width, height)
        if any(box.overlaps(other, 2.0) for other in boxes):
            return False
        boxes.append(box)
    return True


def _marks(construct: Construct, step: int, wanted: int, style: Style) -> list[Mark]:
    return [
        Mark(
            position=tick.position,
            label=tick.label,
            angle=angle_at(tick.position, construct.length, style.start),
        )
        for tick in ruler(construct.length, step, wanted, construct.is_circular)
    ]


def choose_ticks(
    construct: Construct, centre: Point, radius: float, style: Style
) -> tuple[list[Mark], int]:
    """The finest ruler whose labels do not collide, and the step it used."""
    for wanted in TICK_COUNTS:
        step = tick_step(construct.length, wanted)
        marks = _marks(construct, step, wanted, style)
        if _ticks_fit(marks, centre, radius, style):
            return marks, step
    step = tick_step(construct.length, TICK_COUNTS[-1])
    return _marks(construct, step, TICK_COUNTS[-1], style), step


def gc_ring(
    profile: Profile,
    centre: Point,
    baseline: float,
    height: float,
    length: int,
    start: int,
) -> list[str]:
    """The GC track: a curve round a dashed mean, filled between the two.

    The fill is the curve and the mean circle in one path with the even-odd
    rule, which fills exactly the area between them however often the curve
    crosses the mean. The scale is the spread of the windows themselves, so a
    construct whose GC barely moves still shows the shape of what movement
    there is rather than a flat line.
    """
    if profile.highest is None or profile.lowest is None:
        return []
    reach = max((profile.highest - profile.lowest) / 2, 0.01)
    mean = profile.mean or 0.0
    points = []
    for window in profile.windows:
        share = window.gc if window.gc is not None else mean
        offset = (share - mean) / reach * (height / 2)
        points.append(
            polar(
                centre,
                baseline + max(min(offset, height / 2), -height / 2),
                angle_at(window.middle, length, start),
            )
        )
    if not points:
        return []
    curve = polyline(points) + "Z"
    return [
        path(curve + " " + ring(centre, baseline), **{"class": "gc-fill"}),
        path(ring(centre, baseline), **{"class": "gc-axis"}),
        path(curve, **{"class": "gc-line"}),
    ]


def skew_ring(
    profile: Profile,
    centre: Point,
    baseline: float,
    height: float,
    length: int,
    start: int,
) -> list[str]:
    """The GC skew track: one bar per window, coloured by which way it leans.

    Bars rather than a curve, because the useful thing about skew is the sign
    and where it changes, and two colours say that at a glance.
    """
    limit = profile.skew_limit or 1.0
    parts = [path(ring(centre, baseline), **{"class": "gc-axis"})]
    half = profile.step / 2
    for window in profile.windows:
        if window.skew is None:
            continue
        reach = window.skew / limit * (height / 2)
        inner, outer = sorted((baseline, baseline + reach))
        if outer - inner < 0.05:
            continue
        first, last = unwrapped_angles(
            window.middle - half, window.middle + half, length, start
        )
        parts.append(
            path(
                arc_band(centre, inner, outer, first, last),
                **{"class": "skew-up" if window.skew > 0 else "skew-down"},
            )
        )
    return parts


def _features(
    layout: Layout,
    centre: Point,
    bands: list[tuple[float, float]],
    style: Style,
    track: float,
) -> list[str]:
    """Every feature as a band, with the thin line joining separate parts."""
    drawn: list[str] = []
    length = layout.length
    for item in layout.placements:
        inner, outer = bands[item.lane]
        middle = (inner + outer) / 2
        colour = style.palette.colour(item.feature.kind)
        head = min(LONGEST_HEAD, track / max(middle, 1.0) * PER_RADIAN)
        spans = item.spans
        where = (
            f"{item.feature.label} ({item.feature.kind}) "
            f"{location_text(item.feature, length)}"
        )
        for index, span in enumerate(spans):
            first, last = unwrapped_angles(span.start, span.end, length, style.start)
            arrow = "none"
            if item.feature.strand > 0 and index == len(spans) - 1:
                arrow = "end"
            elif item.feature.strand < 0 and index == 0:
                arrow = "start"
            drawn.append(
                path(
                    arc_band(centre, inner, outer, first, last, arrow, head),
                    title=where,
                    fill=colour,
                    **{"class": "feature"},
                )
            )
        for index in range(len(spans) - 1):
            gap_start = spans[index].end
            gap_end = spans[index + 1].start
            if gap_end <= gap_start:
                continue
            first, last = unwrapped_angles(gap_start, gap_end, length, style.start)
            drawn.append(
                path(arc_line(centre, middle, first, last), **{"class": "join"})
            )
    return drawn


def _ruler(construct: Construct, centre: Point, backbone: float, style: Style):
    """The ruler ring, and the step it settled on."""
    label_radius = backbone - TICK - TICK_LABEL_GAP - style.small_size / 2
    marks, step = choose_ticks(construct, centre, label_radius, style)
    if not style.ruler:
        return [], step
    parts = []
    for coordinate in minor_coordinates(construct.length, step):
        angle = angle_at(coordinate - 1, construct.length, style.start)
        parts.append(
            path(
                radial(centre, backbone - TICK / 2, backbone, angle),
                **{"class": "tick-minor"},
            )
        )
    for mark in marks:
        parts.append(
            path(
                radial(centre, backbone - TICK, backbone, mark.angle),
                **{"class": "tick"},
            )
        )
        at = polar(centre, label_radius, mark.angle)
        parts.append(
            text(
                mark.label,
                at.x,
                at.y + style.small_size * 0.35,
                anchor="middle",
                css="tick-label",
            )
        )
    return parts, step


def _labels(placed, centre: Point, radius: float, style: Style) -> list[str]:
    """Every label with the leader line that ties it to its feature."""
    parts: list[str] = []
    for label in placed.labels:
        stub = polar(centre, radius + 3, label.angle)
        elbow = label.x + (-6 if label.side == "right" else 6)
        parts.append(
            path(
                f"M{stub}L{number(elbow)} {number(label.y)}"
                f"L{number(label.x)} {number(label.y)}",
                **{"class": "leader"},
            )
        )
        parts.append(
            text(
                label.text,
                label.x + (2 if label.side == "right" else -2),
                label.y + style.label_size * 0.35,
                anchor="start" if label.side == "right" else "end",
                css="label",
                title=label.full,
                fill=label.colour or None,
            )
        )
    return parts


def draw(
    construct: Construct,
    selection: Selection,
    style: Style,
    profile: Profile | None = None,
) -> Drawing:
    """A circular map of a construct: a standalone SVG and an account of it.

    A construct whose record says it is linear can still be drawn this way -
    there is no reason a cassette cannot be shown as a circle - but nothing
    pretends its ends are joined: the map says so under the name in the
    middle, and so does the summary.
    """
    length = construct.length
    minimum = minimum_width(length)
    layout = assign(selection.kept, length, construct.is_circular, minimum)
    centre = Point(0.0, 0.0)
    lanes = max(layout.lanes, 1)

    requests = (
        label_requests(construct, layout, style) if style.labels != "none" else []
    )
    left_room, right_room = side_room(requests, style)
    gutter = (left_room + right_room) / 2 + (LEADER_GAP if requests else 0.0)
    radius = (style.width - 2 * style.margin) / 2 - gutter
    track, separation = fit_bands(radius, lanes, style)
    backbone = radius - (lanes - 1) * (track + separation) - track
    if backbone < LEAST_RADIUS:
        raise DrawingError(
            f"a map {number(style.width)} wide has no room for {lanes} "
            f"lane{'' if lanes == 1 else 's'} of features beside labels of up "
            f"to {style.label_width} characters; widen it with --width or "
            f"shorten the labels with --label-width"
        )

    body: list[str] = [
        group([path(ring(centre, backbone), **{"class": "backbone"})], "backbone")
    ]
    bands = lane_radii(radius, track, separation, lanes)
    body.append(group(_features(layout, centre, bands, style, track), "features"))

    tick_parts, step = _ruler(construct, centre, backbone, style)
    body.append(group(tick_parts, "ruler"))

    name = style.title or construct.name
    subtitle = f"{length:,} bp"
    reserve = centre_reserve(name, subtitle, style)
    edge = backbone - (
        TICK + TICK_LABEL_GAP + style.small_size if style.ruler else 4.0
    )
    if (style.gc or style.skew) and profile is None:
        raise DrawingError("a GC track needs the sequence, and this input carries none")
    for wanted, make, track_name in (
        (style.gc, gc_ring, "GC"),
        (style.skew, skew_ring, "GC skew"),
    ):
        if not wanted or profile is None:
            continue
        baseline = edge - style.gc_height / 2 - 4
        if baseline - style.gc_height / 2 < reserve:
            raise DrawingError(
                f"there is no room inside the map for a {track_name} track: it "
                f"would reach the {number(reserve)} px kept clear in the middle "
                f"for the name. Widen the map with --width, make the bands "
                f"thinner with --track, or leave a track off"
            )
        body.append(
            group(
                make(profile, centre, baseline, style.gc_height, length, style.start),
                "gc" if track_name == "GC" else "skew",
            )
        )
        edge = baseline - style.gc_height / 2 - 4

    centre_parts = [
        path(ring(centre, reserve), **{"class": "backbone-fill"}),
        text(name, 0, -2, anchor="middle", css="title"),
        text(subtitle, 0, 13, anchor="middle", css="subtitle"),
    ]
    if not construct.is_circular:
        centre_parts.append(
            text("the record calls it linear", 0, 25, anchor="middle", css="subtitle")
        )
    body.append(group(centre_parts, "centre"))

    placed = None
    bounds = Box(-radius, -radius, 2 * radius, 2 * radius)
    if requests:
        placed = place_circular(
            requests,
            centre,
            radius,
            style.label_size,
            gap=LEADER_GAP,
            width=style.label_width,
        )
        body.append(group(_labels(placed, centre, radius, style), "labels"))
        around = placed.bounds()
        if around is not None:
            left = min(bounds.x, around.x)
            top = min(bounds.y, around.y)
            right = max(bounds.right, around.right)
            bottom = max(bounds.bottom, around.bottom)
            bounds = Box(left, top, right - left, bottom - top)

    legend = groups_of(layout, style.palette) if style.legend else []
    if legend:
        parts = []
        for position, (group_name, colour) in enumerate(legend):
            y = bounds.bottom + style.margin + position * LEGEND_ROW
            parts.append(
                rect(bounds.x, y, 10, 10, fill=colour, **{"class": "legend-swatch"})
            )
            parts.append(text(group_name, bounds.x + 14, y + 9, css="legend-label"))
        body.append(group(parts, "legend"))

    extra = (style.margin + len(legend) * LEGEND_ROW) if legend else 0.0
    width = max(style.width, bounds.width + 2 * style.margin)
    height = bounds.height + extra + 2 * style.margin
    shift_x = -bounds.x + (width - bounds.width) / 2
    shift_y = -bounds.y + style.margin
    svg = document(
        [f'<g transform="translate({number(shift_x)} {number(shift_y)})">']
        + body
        + ["</g>"],
        width=width,
        height=height,
        title=f"{name}, {length:,} bp, {construct.topology}",
        description=construct.description,
        label_size=style.label_size,
        small_size=style.small_size,
        title_size=style.title_size,
    )
    return Drawing(
        svg=svg,
        width=width,
        height=height,
        shape="circular",
        construct=construct,
        style=style,
        selection=selection,
        layout=layout,
        placed=placed,
        profile=profile if (style.gc or style.skew) else None,
        tick_step=step,
        track=track,
    )
