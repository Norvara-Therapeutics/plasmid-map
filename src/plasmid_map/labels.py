"""Putting a label beside every feature without putting two in one place.

A map of real annotation has features a few bases apart and features inside
other features, so the labels that belong to them want to be drawn in the same
place. Nothing here drops a label to solve that. On a circular map the labels
are moved along their side of the drawing until none overlaps another and a
leader line is drawn back to the feature; on a linear one they are stacked into
rows above the track. Either way every feature keeps its label, and the summary
says how many had to be moved and how far.

The circular case uses the ordering that minimises how far the labels move in
total: labels are sorted by where they would like to be, runs that collide are
merged into a block, and each block is centred on the average of its members'
wishes. The result is that labels never cross each other - a label whose
feature is higher is drawn higher - which is what makes a leader line
followable.

How wide a label will be is estimated, not measured. There is no font here to
measure with, and embedding one would make a file that is no longer a drawing.
The estimate is a fixed average glyph width, which is what a reader should
know before trusting a map with very long labels: the estimate is used both
for placing labels and for testing that they do not collide, so the two always
agree with each other, and they are both approximations of what a renderer
will do.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from plasmid_map.geometry import Box, Point, polar

#: The average width of a glyph, as a share of the font size, for the
#: sans-serif families the stylesheet asks for.
GLYPH_WIDTH: float = 0.58

#: The pitch between two lines of text, as a share of the font size. This is
#: what labels are moved apart by.
LINE_HEIGHT: float = 1.35

#: How tall the drawn text itself is taken to be, as a share of the font size.
#: It is less than the pitch, which is why two labels one pitch apart do not
#: touch: the difference is the leading between them.
INK_HEIGHT: float = 1.15

#: How many characters a label keeps before it is cut.
LABEL_WIDTH: int = 28

#: What is put at the end of a label that was cut.
ELLIPSIS: str = "\u2026"


class LabelError(ValueError):
    """A set of labels that cannot be placed."""


def shorten(text: str, width: int = LABEL_WIDTH) -> str:
    """A label cut to a width, with an ellipsis where it was cut.

    The whole text is kept in the drawing as a tooltip, so cutting loses
    nothing a reader cannot get back by hovering over the label.
    """
    collapsed = " ".join(text.split())
    if width <= 0 or len(collapsed) <= width:
        return collapsed
    return collapsed[: max(width - 1, 0)] + ELLIPSIS


def estimate_width(text: str, font_size: float) -> float:
    """How wide a piece of text will be drawn, near enough to place it by."""
    return len(text) * font_size * GLYPH_WIDTH


def separate(ideals: Sequence[float], spacing: float) -> list[float]:
    """Positions at least ``spacing`` apart, as near their ideals as possible.

    Labels are sorted by where they want to be; a run that would collide is
    merged into one block and the block is centred on the average of its
    members' wishes; merging repeats until no two blocks collide. This is the
    arrangement with the least total movement that keeps the labels in order,
    and keeping them in order is what makes a leader line followable.
    """
    if not ideals:
        return []
    if spacing <= 0:
        raise LabelError("labels need room between them")
    order = sorted(range(len(ideals)), key=lambda index: (ideals[index], index))
    blocks: list[list[float]] = []
    for index in order:
        blocks.append([ideals[index], 1.0])
        while len(blocks) > 1:
            previous, current = blocks[-2], blocks[-1]
            previous_bottom = previous[0] / previous[1] + previous[1] * spacing / 2
            current_top = current[0] / current[1] - current[1] * spacing / 2
            if current_top < previous_bottom:
                blocks[-2] = [previous[0] + current[0], previous[1] + current[1]]
                blocks.pop()
            else:
                break
    placed: list[float] = []
    for total, count in blocks:
        top = total / count - count * spacing / 2
        placed.extend(top + spacing * (step + 0.5) for step in range(int(count)))
    settled = [0.0] * len(ideals)
    for position, index in zip(placed, order):
        settled[index] = position
    return settled


@dataclass(frozen=True)
class Request:
    """A label that wants to be drawn beside something."""

    text: str
    angle: float
    colour: str = ""
    key: str = ""


@dataclass(frozen=True)
class Label:
    """A label that has been given somewhere to go."""

    text: str
    full: str
    angle: float
    anchor: Point
    x: float
    y: float
    side: str
    box: Box
    colour: str
    key: str
    moved: float

    @property
    def anchored_end(self) -> bool:
        """Whether the text runs leftwards from its anchor point."""
        return self.side == "left"


@dataclass(frozen=True)
class Placed:
    """Every label placed, and what placing them cost."""

    labels: tuple[Label, ...]
    font_size: float

    @property
    def moved(self) -> tuple[Label, ...]:
        """The labels that had to be moved to keep them apart."""
        return tuple(item for item in self.labels if item.moved >= 0.5)

    @property
    def furthest(self) -> float:
        """How far the furthest-moved label moved."""
        return max((item.moved for item in self.labels), default=0.0)

    def bounds(self) -> Box | None:
        """A rectangle holding every label, or None when there are none."""
        if not self.labels:
            return None
        left = min(item.box.x for item in self.labels)
        right = max(item.box.right for item in self.labels)
        top = min(item.box.y for item in self.labels)
        bottom = max(item.box.bottom for item in self.labels)
        return Box(left, top, right - left, bottom - top)

    def collisions(self) -> list[tuple[str, str]]:
        """Any two labels that still overlap, which should be none.

        Kept as a method rather than only as a test so that the drawing can
        say so in its summary if ever it happens.
        """
        found: list[tuple[str, str]] = []
        for first in range(len(self.labels)):
            for second in range(first + 1, len(self.labels)):
                one, two = self.labels[first], self.labels[second]
                if one.box.overlaps(two.box):
                    found.append((one.text, two.text))
        return found


def place_circular(
    requests: Sequence[Request],
    centre: Point,
    radius: float,
    font_size: float,
    gap: float = 16.0,
    width: int = LABEL_WIDTH,
) -> Placed:
    """Labels round a circular map, moved apart along each side.

    A label sits on the side of the drawing its feature is on, at the height
    it would be at if nothing were in the way, and is moved up or down from
    there only as far as it must be. Left-hand labels are drawn ending at
    their column and right-hand ones starting at it, so both read outwards
    from the map.
    """
    pitch = font_size * LINE_HEIGHT
    height = font_size * INK_HEIGHT
    sides: dict[str, list[int]] = {"left": [], "right": []}
    anchors: list[Point] = []
    ideals: list[float] = []
    for index, request in enumerate(requests):
        anchor = polar(centre, radius, request.angle)
        anchors.append(anchor)
        ideals.append(anchor.y)
        turn = request.angle % 360
        sides["right" if turn < 180 else "left"].append(index)
    settled: list[float] = [0.0] * len(requests)
    for side, members in sides.items():
        wanted = [ideals[index] for index in members]
        for index, position in zip(members, separate(wanted, pitch)):
            settled[index] = position
    labels: list[Label] = []
    for side, members in sides.items():
        for index in members:
            request = requests[index]
            text = shorten(request.text, width)
            text_width = estimate_width(text, font_size)
            if side == "right":
                x = centre.x + radius + gap
                box = Box(x, settled[index] - height / 2, text_width, height)
            else:
                x = centre.x - radius - gap
                box = Box(x - text_width, settled[index] - height / 2, text_width, height)
            labels.append(
                Label(
                    text=text,
                    full=" ".join(request.text.split()),
                    angle=request.angle,
                    anchor=anchors[index],
                    x=x,
                    y=settled[index],
                    side=side,
                    box=box,
                    colour=request.colour,
                    key=request.key,
                    moved=abs(settled[index] - ideals[index]),
                )
            )
    labels.sort(key=lambda item: (item.angle, item.text))
    return Placed(labels=tuple(labels), font_size=font_size)


@dataclass(frozen=True)
class Along:
    """A label that wants to sit over a point on a linear track.

    ``target`` is where its leader line should end, which is the top of the
    row its own feature was drawn in rather than the top of the stack: with
    fifteen rows of features, a line stopping at the first of them would not
    say which feature it belonged to.
    """

    text: str
    x: float
    colour: str = ""
    key: str = ""
    target: float = 0.0


@dataclass(frozen=True)
class Stacked:
    """A label on a linear track, with the row it was put in."""

    text: str
    full: str
    x: float
    row: int
    anchor: float
    box: Box
    colour: str
    key: str
    target: float = 0.0


@dataclass(frozen=True)
class Rows:
    """Labels stacked into rows, and how many rows that took."""

    labels: tuple[Stacked, ...]
    rows: int
    font_size: float
    row_height: float
    shortened: int = 0

    def collisions(self) -> list[tuple[str, str]]:
        """Any two labels that still overlap, which should be none."""
        found: list[tuple[str, str]] = []
        for first in range(len(self.labels)):
            for second in range(first + 1, len(self.labels)):
                one, two = self.labels[first], self.labels[second]
                if one.box.overlaps(two.box):
                    found.append((one.text, two.text))
        return found


def stack_linear(
    requests: Sequence[Along],
    top: float,
    font_size: float,
    row_height: float | None = None,
    width: int = LABEL_WIDTH,
    padding: float = 4.0,
    left: float = 0.0,
    right: float | None = None,
) -> Rows:
    """Labels along a linear track, stacked into rows where they collide.

    Each label is centred over its feature and put in the highest row with
    room for it, which keeps a label near what it names rather than moving it
    sideways to somewhere a reader would have to follow a long line to find.
    A label that would run off either end is pulled back inside, since a label
    drawn off the canvas is a label lost.
    """
    height = row_height if row_height is not None else font_size * LINE_HEIGHT
    ink = min(height, font_size * INK_HEIGHT)
    room = None if right is None else max(right - left, 0.0)
    taken: list[list[Box]] = []
    placed: list[Stacked] = []
    shortened = 0
    ordered = sorted(range(len(requests)), key=lambda index: (requests[index].x, index))
    for index in ordered:
        request = requests[index]
        centre_x = request.x
        shown = shorten(request.text, width)
        text_width = estimate_width(shown, font_size)
        if room is not None and text_width > room:
            # Moving a label that is wider than the track cannot bring it
            # inside, so it is cut until it fits. Cutting is visible - the
            # label ends in an ellipsis and keeps its full text as a tooltip -
            # while a label drawn off the canvas is simply lost.
            fits = max(1, int(room / (font_size * GLYPH_WIDTH)))
            shown = shorten(request.text, fits)
            text_width = estimate_width(shown, font_size)
            shortened += 1
        x = centre_x - text_width / 2
        if x < left:
            x = left
        if right is not None and x + text_width > right:
            x = max(left, right - text_width)
        row = 0
        while True:
            if row == len(taken):
                taken.append([])
            box = Box(x, top + row * height + (height - ink) / 2, text_width, ink)
            if not any(box.overlaps(other, padding) for other in taken[row]):
                taken[row].append(box)
                placed.append(
                    Stacked(
                        text=shown,
                        full=" ".join(request.text.split()),
                        x=x + text_width / 2,
                        row=row,
                        anchor=centre_x,
                        box=box,
                        colour=request.colour,
                        key=request.key,
                        target=request.target,
                    )
                )
                break
            row += 1
    placed.sort(key=lambda item: (item.row, item.x))
    return Rows(
        labels=tuple(placed),
        rows=len(taken),
        font_size=font_size,
        row_height=height,
        shortened=shortened,
    )
