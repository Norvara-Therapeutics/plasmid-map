"""Turning base positions into the shapes a drawing is made of.

Angles are measured in degrees from twelve o'clock and increase clockwise,
which is how every plasmid map is read, and the origin of the construct sits
at twelve o'clock unless ``--start`` moves it. Screen coordinates have y
increasing downwards, as SVG does, so the conversion is ``x = sin`` and
``y = -cos`` rather than the other way round; doing it once here is what keeps
every other module free of trigonometry.

Numbers are written to a fixed two decimal places. The same input therefore
draws the same bytes, which is what lets two runs of a map be compared with
``diff`` rather than by eye.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: How many decimal places a coordinate is written to.
PLACES: int = 2

#: Below this many degrees a shape is drawn as a straight-sided wedge, since
#: the arc and the chord differ by less than the width of a drawn line.
FLAT_DEGREES: float = 0.05


class GeometryError(ValueError):
    """A shape that cannot be drawn."""


def number(value: float) -> str:
    """One coordinate, written the one way.

    Negative zero is written as zero: the two are the same point, and letting
    both appear would mean the same map could produce two different files.
    """
    rounded = round(float(value), PLACES)
    if rounded == 0:
        rounded = 0.0
    return f"{rounded:.{PLACES}f}"


@dataclass(frozen=True)
class Point:
    """A point on the canvas."""

    x: float
    y: float

    def __str__(self) -> str:
        return f"{number(self.x)},{number(self.y)}"


@dataclass(frozen=True)
class Box:
    """A rectangle, used for deciding whether two labels collide."""

    x: float
    y: float
    width: float
    height: float

    @property
    def right(self) -> float:
        """The far edge."""
        return self.x + self.width

    @property
    def bottom(self) -> float:
        """The lower edge."""
        return self.y + self.height

    def overlaps(self, other: Box, padding: float = 0.0) -> bool:
        """Whether two rectangles share any area, with room left between."""
        return (
            self.x - padding < other.right
            and other.x - padding < self.right
            and self.y - padding < other.bottom
            and other.y - padding < self.bottom
        )

    def grown(self, by: float) -> Box:
        """The same rectangle with a margin round it."""
        return Box(self.x - by, self.y - by, self.width + 2 * by, self.height + 2 * by)


def angle_at(position: float, length: int, start: int = 0) -> float:
    """The angle, in degrees from twelve o'clock, of a base position.

    ``start`` is the base put at twelve o'clock, which is how a map is rotated
    so that a feature of interest sits at the top.
    """
    if length <= 0:
        raise GeometryError("a construct with no length has no angles")
    return ((position - start) % length) * 360.0 / length


def unwrapped_angles(
    first: float, last: float, length: int, start: int = 0
) -> tuple[float, float]:
    """The angles of a span's ends, with the second kept after the first.

    A span that crosses twelve o'clock would otherwise end at a smaller angle
    than it began at, and be drawn the long way round the circle.
    """
    a = angle_at(first, length, start)
    b = a + (last - first) * 360.0 / length
    return a, b


def polar(centre: Point, radius: float, angle: float) -> Point:
    """The point at a radius and an angle from twelve o'clock, clockwise."""
    radians = math.radians(angle)
    return Point(
        centre.x + radius * math.sin(radians),
        centre.y - radius * math.cos(radians),
    )


def _arc(
    centre: Point, radius: float, start: float, end: float, forwards: bool
) -> str:
    """One elliptical-arc command from one angle to another."""
    target = polar(centre, radius, end)
    large = 1 if abs(end - start) > 180 else 0
    sweep = 1 if forwards else 0
    return (
        f"A{number(radius)} {number(radius)} 0 {large} {sweep} "
        f"{number(target.x)} {number(target.y)}"
    )


def ring(centre: Point, radius: float) -> str:
    """A whole circle as a path, so that it strokes like every other shape."""
    top = polar(centre, radius, 0)
    bottom = polar(centre, radius, 180)
    return (
        f"M{number(top.x)} {number(top.y)}"
        f"A{number(radius)} {number(radius)} 0 1 1 "
        f"{number(bottom.x)} {number(bottom.y)}"
        f"A{number(radius)} {number(radius)} 0 1 1 "
        f"{number(top.x)} {number(top.y)}Z"
    )


def arc_band(
    centre: Point,
    inner: float,
    outer: float,
    start: float,
    end: float,
    arrow: str = "none",
    arrow_degrees: float = 0.0,
) -> str:
    """A band between two radii, between two angles, with an optional point.

    ``arrow`` is ``start``, ``end`` or ``none``, and draws the band's leading
    edge as a point rather than a straight cut, which is how a feature's
    strand is shown. The point never takes more than the band's own width, so
    a short feature becomes a triangle rather than an arrowhead pointing back
    past its own start.
    """
    if outer <= inner:
        raise GeometryError("a band needs an outer radius larger than its inner")
    if end < start:
        start, end = end, start
    span = end - start
    if span >= 360:
        return ring(centre, outer) + " " + ring(centre, inner)
    middle = (inner + outer) / 2
    head = max(0.0, min(arrow_degrees, span))
    if arrow in ("start", "end") and head >= span > 0:
        # The arrowhead is as long as the feature, so the feature is the
        # arrowhead: a triangle, with no arc left to draw either side of it.
        back, point = (end, start) if arrow == "start" else (start, end)
        return (
            f"M{polar(centre, outer, back)}"
            f"L{polar(centre, middle, point)}"
            f"L{polar(centre, inner, back)}Z"
        )
    if arrow == "end" and head > 0:
        outer_end = end - head
        tip = polar(centre, middle, end)
        parts = [
            f"M{polar(centre, outer, start)}",
            _arc(centre, outer, start, outer_end, True),
            f"L{tip}",
            f"L{polar(centre, inner, outer_end)}",
            _arc(centre, inner, outer_end, start, False),
            "Z",
        ]
    elif arrow == "start" and head > 0:
        outer_start = start + head
        tip = polar(centre, middle, start)
        parts = [
            f"M{tip}",
            f"L{polar(centre, outer, outer_start)}",
            _arc(centre, outer, outer_start, end, True),
            f"L{polar(centre, inner, end)}",
            _arc(centre, inner, end, outer_start, False),
            "Z",
        ]
    else:
        parts = [
            f"M{polar(centre, outer, start)}",
            _arc(centre, outer, start, end, True),
            f"L{polar(centre, inner, end)}",
            _arc(centre, inner, end, start, False),
            "Z",
        ]
    return "".join(parts)


def arc_line(centre: Point, radius: float, start: float, end: float) -> str:
    """A single arc with no thickness, for a tick or a joining line."""
    if end < start:
        start, end = end, start
    if end - start >= 360:
        return ring(centre, radius)
    return f"M{polar(centre, radius, start)}" + _arc(centre, radius, start, end, True)


def radial(centre: Point, inner: float, outer: float, angle: float) -> str:
    """A straight line out from the centre, for a tick mark."""
    return f"M{polar(centre, inner, angle)}L{polar(centre, outer, angle)}"


def chord_length(radius: float, degrees: float) -> float:
    """How long a span is along the page, for deciding whether a label fits."""
    return 2 * radius * math.sin(math.radians(min(degrees, 180.0)) / 2)


def polyline(points: list[Point]) -> str:
    """A path through a list of points, for a profile track."""
    if not points:
        return ""
    first, *rest = points
    return f"M{first}" + "".join(f"L{point}" for point in rest)
