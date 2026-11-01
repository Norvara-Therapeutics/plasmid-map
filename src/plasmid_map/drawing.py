"""What the two renderers share: the options they take and what they report.

A drawing is not only a file. It is also an account of what was drawn and what
was not: which features were left out and why, which were drawn wider than
they are, how many labels had to be moved, which feature kinds had no colour
of their own. A map that quietly dropped a feature would look exactly like a
map of a construct that did not have one, so nothing is dropped quietly.

Selection happens here rather than in each renderer, so a circular map and a
linear one of the same file show the same features.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

from plasmid_map.composition import Profile
from plasmid_map.feature import Construct, Feature
from plasmid_map.genbank import WHOLE_RECORD_KINDS
from plasmid_map.labels import Placed, Rows
from plasmid_map.lanes import Layout
from plasmid_map.palette import Palette, is_default

#: How a label may be placed.
LABEL_MODES: tuple[str, ...] = ("outside", "none")

#: The width of the canvas, in user units, before labels widen it.
WIDTH: float = 760.0


class DrawingError(ValueError):
    """A drawing that cannot be made as asked."""


@dataclass(frozen=True)
class Style:
    """Everything about how a map looks, with the defaults it is drawn at."""

    palette: Palette = field(default_factory=Palette)
    width: float = WIDTH
    track: float = 15.0
    separation: float = 4.0
    margin: float = 18.0
    label_size: float = 11.0
    small_size: float = 9.0
    title_size: float = 15.0
    labels: str = "outside"
    label_width: int = 28
    ruler: bool = True
    gc: bool = False
    skew: bool = False
    gc_height: float = 30.0
    legend: bool = False
    title: str = ""
    start: int = 0

    def __post_init__(self) -> None:
        if self.labels not in LABEL_MODES:
            raise DrawingError(
                f"not a way to place labels: {self.labels!r}; "
                f"it is one of {', '.join(LABEL_MODES)}"
            )
        if self.width <= 0:
            raise DrawingError(f"a canvas needs a width, not {self.width}")
        if self.track <= 0:
            raise DrawingError(f"a feature band needs a thickness, not {self.track}")


@dataclass(frozen=True)
class Selection:
    """The features that will be drawn, and an account of the rest."""

    kept: tuple[Feature, ...]
    hidden: tuple[tuple[str, int], ...]
    whole_record: tuple[Feature, ...]
    short: tuple[Feature, ...]
    unknown_kinds: tuple[str, ...]

    @property
    def dropped(self) -> int:
        """How many features are not being drawn."""
        return (
            sum(count for _kind, count in self.hidden)
            + len(self.whole_record)
            + len(self.short)
        )


def select(
    construct: Construct,
    kinds: Sequence[str] = (),
    hide: Sequence[str] = (),
    minimum_length: int = 0,
    whole_record: bool = False,
) -> Selection:
    """Which features to draw, and why the others are not drawn.

    By default the feature that describes the whole record - ``source`` in a
    GenBank file - is left out, since a band round the entire construct says
    nothing and hides a lane. Asking for it by kind brings it back.
    """
    wanted = set(kinds)
    unwanted = set(hide)
    kept: list[Feature] = []
    skipped: dict[str, int] = {}
    whole: list[Feature] = []
    short: list[Feature] = []
    for feature in construct.features:
        if feature.kind in unwanted or (wanted and feature.kind not in wanted):
            skipped[feature.kind] = skipped.get(feature.kind, 0) + 1
            continue
        if feature.kind in WHOLE_RECORD_KINDS and not (
            whole_record or feature.kind in wanted
        ):
            whole.append(feature)
            continue
        if minimum_length and feature.length < minimum_length:
            short.append(feature)
            continue
        kept.append(feature)
    unknown = sorted({feature.kind for feature in kept if is_default(feature.kind)})
    return Selection(
        kept=tuple(kept),
        hidden=tuple(sorted(skipped.items())),
        whole_record=tuple(whole),
        short=tuple(short),
        unknown_kinds=tuple(unknown),
    )


def groups_of(layout: Layout, palette: Palette) -> list[tuple[str, str]]:
    """What a legend should say: one row per colour actually drawn.

    Groups rather than kinds, because a legend naming eleven kinds of binding
    site teaches a reader nothing a legend naming 'binding' does not - except
    for a kind given a colour of its own with ``--colour``, which gets its own
    row under its own name. Labelling an overridden kind with its group would
    put one swatch against two colours and tell the reader the wrong one.
    """
    from plasmid_map.palette import group_of

    seen: dict[tuple[str, str], None] = {}
    for item in layout.placements:
        kind = item.feature.kind
        label = kind if kind in palette.overrides else group_of(kind)
        seen.setdefault((label, palette.colour(kind)), None)
    return sorted(seen)


@dataclass(frozen=True)
class Drawing:
    """A finished map and the account of how it was drawn."""

    svg: str
    width: float
    height: float
    shape: str
    construct: Construct
    style: Style
    selection: Selection
    layout: Layout
    placed: Placed | Rows | None
    profile: Profile | None = None
    tick_step: int = 0
    track: float = 0.0
    split: tuple[Feature, ...] = ()
    added: int = 0

    @property
    def thinned(self) -> bool:
        """Whether the feature bands were made thinner than they were asked.

        They are, when a construct needs more lanes than the map has room for
        at the full thickness. The alternative is lanes reaching the centre,
        which leaves nowhere for the ruler and the name.
        """
        return bool(self.track) and self.track < self.style.track - 0.01

    @property
    def drawn(self) -> int:
        """How many features are in the drawing."""
        return len(self.layout.placements)

    @property
    def crossing(self) -> tuple[Feature, ...]:
        """The drawn features that run past the origin."""
        return tuple(
            item.feature
            for item in self.layout.placements
            if item.feature.crosses_origin(self.construct.length)
        )

    @property
    def partial(self) -> tuple[Feature, ...]:
        """The drawn features with an end that runs off the record."""
        return tuple(
            item.feature for item in self.layout.placements if item.feature.is_partial
        )

    @property
    def label_collisions(self) -> list[tuple[str, str]]:
        """Any two labels still overlapping, which should always be none."""
        if self.placed is None:
            return []
        return self.placed.collisions()

    def groups(self) -> list[tuple[str, str]]:
        """The feature groups present and their colours, for a legend."""
        return groups_of(self.layout, self.style.palette)
