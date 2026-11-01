"""Where the tick marks go and what they say.

Ticks are placed on a 1-2-5 scale - 100, 200, 500, 1,000 - because those are
the intervals a reader adds up without thinking. The step chosen is the
coarsest that still gives at least as many ticks as were asked for, so a 2,686
base plasmid gets a tick every 200 bases rather than every 268.

Tick positions count from 1, as a person reading a sequence does, while
everything inside the program counts from 0. The conversion happens here and
nowhere else.
"""

from __future__ import annotations

from dataclasses import dataclass

#: The mantissas a step may take, in the order they are tried.
STEPS: tuple[int, ...] = (1, 2, 5)

#: How many ticks a ruler aims for before the scale is chosen.
WANTED: int = 10

#: How many minor divisions sit between two major ticks.
MINOR_DIVISIONS: int = 5


class RulerError(ValueError):
    """A ruler that cannot be drawn."""


def tick_step(length: int, wanted: int = WANTED) -> int:
    """The spacing between major ticks for a construct of this length.

    The coarsest 1-2-5 step that still gives ``wanted`` ticks or more, and
    never zero, so a very short construct is ticked every base rather than not
    at all.
    """
    if length <= 0:
        raise RulerError("a construct with no length has no ruler")
    if wanted <= 0:
        raise RulerError("a ruler needs at least one tick")
    rough = length / wanted
    if rough <= 1:
        return 1
    magnitude = 1
    while magnitude * 10 <= rough:
        magnitude *= 10
    step = magnitude
    for mantissa in STEPS:
        if magnitude * mantissa <= rough:
            step = magnitude * mantissa
    return step


def coordinates(
    length: int, step: int | None = None, wanted: int = WANTED, circular: bool = True
) -> list[int]:
    """The 1-based coordinates the major ticks mark.

    Base 1 is always marked, since it is the origin and a map has to say where
    the numbering starts. A tick that would land on the last base of a
    circular construct is left out, because it is the same place as base 1 and
    would put two labels on one mark.
    """
    spacing = tick_step(length, wanted) if step is None else step
    if spacing <= 0:
        raise RulerError(f"a tick step has to move along, not {spacing}")
    marked = [1]
    position = spacing
    while position <= length:
        if position != 1 and not (circular and position == length):
            marked.append(position)
        position += spacing
    return marked


def minor_coordinates(
    length: int, step: int, divisions: int = MINOR_DIVISIONS
) -> list[int]:
    """The 1-based coordinates of the small ticks between the major ones."""
    if divisions <= 1:
        return []
    small = max(1, step // divisions)
    return [
        position
        for position in range(small, length + 1, small)
        if position % step
    ]


@dataclass(frozen=True)
class Tick:
    """One tick: the 0-based position it is drawn at, and what it says."""

    position: int
    coordinate: int
    label: str


def tick_label(coordinate: int, step: int) -> str:
    """What is written at a tick.

    Thousands become kilobases once the step is a round number of them, which
    keeps a long construct's ruler from being a row of five-digit numbers.
    Base 1 keeps its number either way, since it is the one tick whose
    position is a convention rather than a measurement.
    """
    if coordinate == 1:
        return "1"
    if step >= 1000 and step % 1000 == 0 and coordinate % 1000 == 0:
        return f"{coordinate // 1000:,} kb"
    return f"{coordinate:,}"


def ruler(
    length: int, step: int | None = None, wanted: int = WANTED, circular: bool = True
) -> list[Tick]:
    """Every major tick, where it is drawn and what it says."""
    spacing = tick_step(length, wanted) if step is None else step
    return [
        Tick(
            position=coordinate - 1,
            coordinate=coordinate,
            label=tick_label(coordinate, spacing),
        )
        for coordinate in coordinates(length, spacing, wanted, circular)
    ]
