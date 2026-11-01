"""What is drawn: the spans a feature occupies and the construct beneath it.

Positions are held 0-based and half-open, whatever they were read from.
``start`` is the first base and ``end`` the base after the last, so a feature
written ``86..1276`` in a GenBank file is ``Span(85, 1276)`` here and is 1,191
bases long. GenBank counts from 1 with both ends included and is converted on
the way in, once, in the reader, so nothing downstream sees two conventions.

A feature is one or more parts. Separate parts are separate pieces of one
feature - the exons of a spliced coding sequence - drawn as blocks joined by a
thin line. A part whose ``end`` runs past the length of a circular construct is
a single piece that crosses the origin, and is drawn as one arc. Holding it
that way rather than as two parts is what keeps its length right: phiX174's A
gene is written ``join(3981..5386,1..136)`` and is 1,542 bases long, not the
5,386 that subtracting its first coordinate from its last would give, and not
the two features that drawing each half separately would show.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

#: The topologies a construct may have.
TOPOLOGIES: tuple[str, ...] = ("circular", "linear")

#: The strands, as a sign. 0 is a feature with no strand of its own.
STRANDS: tuple[int, ...] = (1, -1, 0)

#: How a strand is written in a table and read back from one.
STRAND_SYMBOLS: dict[int, str] = {1: "+", -1: "-", 0: "."}


class FeatureError(ValueError):
    """A feature or construct that cannot be held as it stands."""


def strand_of(symbol: str) -> int:
    """The sign of a strand written as a symbol.

    ``+`` and ``-`` are the signs; anything a format uses for 'no strand of
    its own' - ``.``, ``?``, an empty field - is 0.
    """
    text = symbol.strip()
    if text in ("+", "1", "+1"):
        return 1
    if text in ("-", "-1"):
        return -1
    if text in ("", ".", "?", "0"):
        return 0
    raise FeatureError(f"not a strand: {symbol!r}")


@dataclass(frozen=True, order=True)
class Span:
    """A half-open run of bases, which may run past the origin.

    ``end`` is allowed to exceed the length of a circular construct: a feature
    crossing the origin of a 5,386 base plasmid is one span from 3,980 to
    5,522, not two spans whose lengths a later sum would have to know not to
    add up differently. ``end`` equal to ``start`` is a site between two bases,
    which GenBank writes ``102^103``; it has no length and is drawn as a mark.
    """

    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise FeatureError(f"a span cannot start before base 1: {self.start}")
        if self.end < self.start:
            raise FeatureError(f"a span cannot end before it starts: {self!r}")

    @property
    def length(self) -> int:
        """How many bases the span holds; zero for a site between two."""
        return self.end - self.start

    def shifted(self, by: int) -> Span:
        """The same span moved along the sequence."""
        return Span(self.start + by, self.end + by)

    def widened(self, to: int) -> Span:
        """The span grown about its middle until it is at least this long.

        Drawing uses this and so does lane assignment, which is the point: a
        one base feature drawn as wide as twenty has to be kept apart from its
        neighbours as though it were twenty, or two marks a few bases apart
        would be drawn on top of each other.
        """
        if self.length >= to:
            return self
        middle = (self.start + self.end) / 2
        start = int(round(middle - to / 2))
        return Span(max(start, 0), max(start, 0) + to)


def overlapping(
    first: Span, second: Span, length: int, circular: bool, padding: int = 0
) -> bool:
    """Whether two spans share a base, counting round the origin.

    On a circular construct a span near the end and a span near the start may
    be neighbours, so the second is also tried shifted a whole turn each way.
    The shift is done in arithmetic rather than by building a span, since a
    span shifted back a turn starts before the first base and there is no such
    position. ``padding`` grows the first span on both sides, which is how a
    lane keeps clear space between what it holds.

    Spans that merely touch do not overlap, as half-open spans do not.
    """
    low = first.start - padding
    high = first.end + padding
    offsets = (-length, 0, length) if circular else (0,)
    for offset in offsets:
        if low < second.end + offset and second.start + offset < high:
            return True
    return False


@dataclass(frozen=True)
class Feature:
    """One annotated thing on a construct.

    ``kind`` is the GenBank feature key - CDS, gene, promoter, rep_origin -
    which is what decides the colour and the shape it is drawn with. ``name``
    is what the label says.
    """

    kind: str
    name: str = ""
    parts: tuple[Span, ...] = ()
    strand: int = 0
    qualifiers: tuple[tuple[str, str], ...] = ()
    partial: tuple[bool, bool] = (False, False)

    def __post_init__(self) -> None:
        if not self.kind:
            raise FeatureError("a feature needs a kind")
        if not self.parts:
            raise FeatureError(f"{self.kind} has no position")
        if self.strand not in STRANDS:
            raise FeatureError(f"not a strand: {self.strand!r}")
        if list(self.parts) != sorted(self.parts):
            object.__setattr__(self, "parts", tuple(sorted(self.parts)))

    @property
    def start(self) -> int:
        """The first base of the first part."""
        return self.parts[0].start

    @property
    def end(self) -> int:
        """The base after the last base of the last part."""
        return max(part.end for part in self.parts)

    @property
    def extent(self) -> Span:
        """First base to last, parts and the gaps between them together."""
        return Span(self.start, self.end)

    @property
    def length(self) -> int:
        """How many bases the feature holds, gaps between parts excluded.

        This is the sum of the parts, never the distance from the first
        coordinate to the last, which for a spliced or origin-crossing feature
        are different numbers.
        """
        return sum(part.length for part in self.parts)

    @property
    def label(self) -> str:
        """What to write beside it: its name, or its kind when it has none."""
        return self.name or self.kind

    def qualifier(self, key: str) -> str | None:
        """One qualifier's value, or None when the feature has no such key."""
        for name, value in self.qualifiers:
            if name == key:
                return value
        return None

    def crosses_origin(self, length: int) -> bool:
        """Whether any part of it runs past the end of the construct."""
        return any(part.end > length for part in self.parts)

    @property
    def is_partial(self) -> bool:
        """Whether either end of it is known to run past what was sequenced.

        GenBank writes those ends ``<1`` and ``>5386``. The flag is kept and
        reported rather than quietly rounded to the end of the record, since a
        feature that continues off the end of a construct and a feature that
        stops exactly at it are different claims.
        """
        return any(self.partial)

    def with_parts(self, parts: Sequence[Span]) -> Feature:
        """The same feature at different positions."""
        return Feature(
            kind=self.kind,
            name=self.name,
            parts=tuple(parts),
            strand=self.strand,
            qualifiers=self.qualifiers,
            partial=self.partial,
        )


def joined(parts: Sequence[Span], length: int, circular: bool) -> tuple[Span, ...]:
    """Parts with any that meet at the origin made into one.

    ``join(3981..5386,1..136)`` on a circular molecule is one piece that
    crosses the origin, while the same join on a linear molecule is two pieces
    with a gap between them. Only the first reading is merged, and only when
    the parts actually meet: a join of ``4000..5386`` and ``50..100`` leaves a
    gap and stays two parts.
    """
    ordered = list(parts)
    if not circular or len(ordered) < 2:
        return tuple(sorted(ordered))
    merged: list[Span] = []
    tail = [part for part in ordered if part.end == length]
    head = [part for part in ordered if part.start == 0]
    if len(tail) == 1 and len(head) == 1 and tail[0] is not head[0]:
        crossing = Span(tail[0].start, length + head[0].length)
        merged = [part for part in ordered if part not in (tail[0], head[0])]
        merged.append(crossing)
        return tuple(sorted(merged))
    return tuple(sorted(ordered))


def within(feature: Feature, length: int, circular: bool) -> bool:
    """Whether every part of a feature lies on a construct of this length.

    A part of a circular construct may run past the end, by no more than one
    whole turn. A part of a linear one may not run past the end at all.
    """
    for part in feature.parts:
        if circular and part.length == 0 and part.start == length:
            # ``5386^1`` is the site between the last base and the first,
            # which on a circular molecule is a real place: the origin.
            continue
        if part.start >= length and not (part.start == 0 and length == 0):
            return False
        if circular:
            if part.end > part.start + length:
                return False
        elif part.end > length:
            return False
    return True


@dataclass(frozen=True)
class Construct:
    """A sequence and what is annotated on it.

    The sequence is optional: a feature table says where things are without
    saying what the bases are, and a map can be drawn from either. What cannot
    be drawn from a table alone is the GC track, and asking for one without a
    sequence is refused rather than drawn empty.
    """

    name: str
    length: int
    topology: str = "circular"
    features: tuple[Feature, ...] = ()
    sequence: str | None = None
    description: str = ""
    accession: str = ""
    molecule: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        if self.length <= 0:
            raise FeatureError(f"a construct needs a length, not {self.length}")
        if self.topology not in TOPOLOGIES:
            raise FeatureError(
                f"not a topology: {self.topology!r}; "
                f"it is one of {', '.join(TOPOLOGIES)}"
            )
        if self.sequence is not None and len(self.sequence) != self.length:
            raise FeatureError(
                f"{self.name} declares {self.length:,} bases and carries "
                f"{len(self.sequence):,}"
            )
        for feature in self.features:
            if not within(feature, self.length, self.is_circular):
                raise FeatureError(
                    f"{feature.label} lies at {feature.start + 1}..{feature.end} "
                    f"on a {self.topology} construct of {self.length:,} bases"
                )

    @property
    def is_circular(self) -> bool:
        """Whether the ends of the sequence are joined."""
        return self.topology == "circular"

    def kinds(self) -> list[tuple[str, int]]:
        """The feature kinds present, most common first then alphabetical."""
        counted: dict[str, int] = {}
        for feature in self.features:
            counted[feature.kind] = counted.get(feature.kind, 0) + 1
        return sorted(counted.items(), key=lambda pair: (-pair[1], pair[0]))

    def crossing(self) -> tuple[Feature, ...]:
        """The features that run past the origin."""
        return tuple(
            feature
            for feature in self.features
            if feature.crosses_origin(self.length)
        )

    def with_features(self, features: Iterable[Feature]) -> Construct:
        """The same construct carrying a different set of features."""
        return Construct(
            name=self.name,
            length=self.length,
            topology=self.topology,
            features=tuple(features),
            sequence=self.sequence,
            description=self.description,
            accession=self.accession,
            molecule=self.molecule,
            source=self.source,
        )


def location_text(feature: Feature, length: int) -> str:
    """Where a feature lies, written the way a GenBank file writes it.

    One place produces this string and everything that shows a position uses
    it: the ``location`` column of every table, and the tooltip on every shape
    in a drawing. A site between two bases comes out ``102^103`` rather than
    as a pair of coordinates that would have to be read as a range, and a
    feature crossing the origin comes out as the join it was written as.
    """
    pieces: list[str] = []
    for part in feature.parts:
        if part.end > length:
            pieces.append(f"{part.start + 1}..{length}")
            pieces.append(f"1..{part.end - length}")
        elif part.length == 0:
            after = part.start + 1
            pieces.append(f"{part.start}^{1 if after > length else after}")
        elif part.length == 1:
            pieces.append(f"{part.start + 1}")
        else:
            pieces.append(f"{part.start + 1}..{part.end}")
    # A complement location is written with its markers in the order the
    # record wrote them, which is the reverse of the order they are held in:
    # ``complement(<1..300)`` parses to a feature whose *last* end is partial.
    first_partial, last_partial = feature.partial
    if feature.strand < 0:
        first_partial, last_partial = last_partial, first_partial
    if first_partial and pieces:
        pieces[0] = "<" + pieces[0]
    if last_partial and pieces:
        head, _, tail = pieces[-1].rpartition("..")
        pieces[-1] = f"{head}..>{tail}" if head else ">" + pieces[-1]
    written = pieces[0] if len(pieces) == 1 else f"join({','.join(pieces)})"
    return f"complement({written})" if feature.strand < 0 else written


def sorted_features(features: Iterable[Feature]) -> list[Feature]:
    """Features in position order, longest first where they start together.

    Sorting longest first is what puts an enclosing feature outside the ones
    it contains when lanes are assigned, which is how a reader expects to see
    a gene and the coding sequence inside it.
    """
    return sorted(
        features, key=lambda item: (item.start, -item.length, item.kind, item.name)
    )


def minimum_width(length: int, share: float = 0.002) -> int:
    """How wide a feature is drawn at the least, in bases.

    A single base of a 5,000 base plasmid is seven hundredths of a degree and
    would not appear at all, so everything shorter than this is drawn this
    wide. The summary says how many features that happened to, because a mark
    drawn wider than the thing it marks is a thing a reader should be told
    about rather than left to assume.
    """
    return max(1, int(round(length * share)))


def drawn_spans(
    feature: Feature, length: int, minimum: int
) -> tuple[Span, ...]:
    """The spans a feature is actually drawn as, nothing narrower than minimum.

    A feature whose parts are each too short is widened part by part; one
    whose parts are short but which is long overall - a spliced coding
    sequence of many small exons - keeps each part widened as well, because
    each part is drawn.
    """
    return tuple(part.widened(minimum) for part in feature.parts)
