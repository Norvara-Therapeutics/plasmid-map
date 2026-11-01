"""Keeping features that overlap from being drawn on top of each other.

Real annotation overlaps constantly: a gene and the coding sequence inside it,
a promoter inside the region it drives, three phage genes reading out of the
same stretch of DNA in different frames. Drawing them all on one ring would
hide whichever was drawn first, so each feature is given a lane, and lanes are
drawn as concentric rings on a circular map and as stacked rows on a linear
one.

Lane 0 is the outermost ring. Features are placed longest first where they
start together, so an enclosing feature takes the outer lane and what it
contains sits inside it, which is the arrangement a reader expects.

Overlap is measured on the span a feature is actually drawn as, not on the
span it occupies. A one base site drawn as wide as twenty has to be kept
twenty bases clear of its neighbours, or two marks a few bases apart would be
drawn over each other however carefully their true positions were kept.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from plasmid_map.feature import (
    Feature,
    Span,
    drawn_spans,
    overlapping,
    sorted_features,
)

#: How many bases of clear space are left between two features in one lane,
#: as a share of the construct. Without it, two features that merely touch
#: would be drawn as one block.
GAP_SHARE: float = 0.0015


@dataclass(frozen=True)
class Placement:
    """One feature, the lane it was given and the spans it is drawn as."""

    feature: Feature
    lane: int
    spans: tuple[Span, ...]

    @property
    def extent(self) -> Span:
        """First drawn base to last, the gaps between parts included."""
        return Span(self.spans[0].start, max(span.end for span in self.spans))

    @property
    def widened(self) -> bool:
        """Whether the feature is drawn wider than it really is."""
        drawn = sum(span.length for span in self.spans)
        return drawn > self.feature.length


@dataclass(frozen=True)
class Layout:
    """Every feature placed, and how many lanes that took."""

    placements: tuple[Placement, ...]
    lanes: int
    minimum: int
    length: int
    circular: bool

    @property
    def widened(self) -> tuple[Placement, ...]:
        """The features drawn wider than they are."""
        return tuple(item for item in self.placements if item.widened)

    def in_lane(self, lane: int) -> tuple[Placement, ...]:
        """Everything placed in one lane, in position order."""
        return tuple(item for item in self.placements if item.lane == lane)


def gap_for(length: int) -> int:
    """How many bases of clear space a lane keeps between its features."""
    return max(1, int(round(length * GAP_SHARE)))


def assign(
    features: Iterable[Feature],
    length: int,
    circular: bool,
    minimum: int,
    gap: int | None = None,
) -> Layout:
    """Give every feature a lane in which nothing else overlaps it.

    Each feature takes the outermost lane with room for it, which keeps the
    number of lanes to the smallest the arrangement allows and keeps a map
    from growing a ring for a feature that would have fitted beside another.
    """
    ordered = sorted_features(features)
    spacing = gap_for(length) if gap is None else gap
    occupied: list[list[Span]] = []
    placements: list[Placement] = []
    for feature in ordered:
        spans = drawn_spans(feature, length, minimum)
        extent = Span(spans[0].start, max(span.end for span in spans))
        lane = 0
        while True:
            if lane == len(occupied):
                occupied.append([])
            if not any(
                overlapping(extent, taken, length, circular, spacing)
                for taken in occupied[lane]
            ):
                occupied[lane].append(extent)
                break
            lane += 1
        placements.append(Placement(feature=feature, lane=lane, spans=spans))
    return Layout(
        placements=tuple(placements),
        lanes=len(occupied),
        minimum=minimum,
        length=length,
        circular=circular,
    )


def lane_radii(
    outer: float, thickness: float, separation: float, lanes: int
) -> list[tuple[float, float]]:
    """The inner and outer radius of each lane, outermost lane first."""
    radii: list[tuple[float, float]] = []
    edge = outer
    for _lane in range(max(lanes, 1)):
        radii.append((edge - thickness, edge))
        edge -= thickness + separation
    return radii


def lane_rows(
    top: float, thickness: float, separation: float, lanes: int
) -> list[tuple[float, float]]:
    """The upper and lower edge of each row, first lane at the top."""
    rows: list[tuple[float, float]] = []
    edge = top
    for _lane in range(max(lanes, 1)):
        rows.append((edge, edge + thickness))
        edge += thickness + separation
    return rows


def deepest(placements: Sequence[Placement]) -> int:
    """How many lanes a set of placements actually uses."""
    return max((item.lane for item in placements), default=-1) + 1
