"""Reading GenBank flat files, for the three things a map needs from one.

A map needs how long the construct is, whether its ends are joined, and where
its features lie. Everything else in a GenBank record - the references, the
comment block, the taxonomy - is skipped rather than modelled, because this
program never writes a GenBank file and so has nothing it could lose.

The LOCUS line is read left to right by column, not right to left by counting
back from the end: a record whose line stops after the molecule type has no
division and no date, and reading backwards would take its topology from
whatever happened to be last. Topology comes from that line and decides
whether a join at the origin is one feature crossing it or two pieces of a
spliced one, so getting it from the right field matters.

What is refused rather than guessed: a location in another record
(``J00194.1:100..202``), which is not on this construct and cannot be drawn on
it, and ``one-of(...)``, which names several positions where a drawing can
show only one.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterator, Sequence

from plasmid_map.feature import (
    Construct,
    Feature,
    FeatureError,
    Span,
    joined,
)

#: Where a feature key starts and where its location does, in a FEATURES
#: block. Both are fixed by the format; qualifiers are indented to the second.
KEY_COLUMN: int = 5
LOCATION_COLUMN: int = 21

#: Which qualifier gives a feature its label, most preferred first. A gene
#: symbol is shorter and better known than the product description that
#: usually sits beside it, which is why it comes first.
NAME_QUALIFIERS: tuple[str, ...] = (
    "gene",
    "label",
    "locus_tag",
    "standard_name",
    "product",
    "bound_moiety",
    "mobile_element_type",
    "rpt_family",
    "organism",
    "note",
)

#: Feature keys that describe the record rather than a thing on it.
WHOLE_RECORD_KINDS: frozenset[str] = frozenset({"source"})


class GenBankError(ValueError):
    """A GenBank record that cannot be read."""


@dataclass(frozen=True)
class Location:
    """Where a feature lies, as the location string said it.

    ``parts`` are 0-based and half-open. ``partial`` says whether the record
    wrote ``<`` before the first position or ``>`` before the last, meaning
    the feature runs past what was sequenced.
    """

    parts: tuple[Span, ...]
    strand: int
    partial: tuple[bool, bool]


def _position(text: str) -> tuple[int, bool]:
    """One coordinate, counting from 1, and whether it was marked partial."""
    stripped = text.strip()
    partial = stripped.startswith("<") or stripped.startswith(">")
    digits = stripped.lstrip("<>")
    if not digits.isdigit():
        raise GenBankError(f"not a position: {text!r}")
    return int(digits), partial


def _simple(text: str) -> Location:
    """A single position, a range, or a site between two bases."""
    if ":" in text:
        raise GenBankError(
            f"{text!r} points at another record, which is not on this construct"
        )
    if ".." in text:
        left, _, right = text.partition("..")
        start, left_partial = _position(left)
        end, right_partial = _position(right)
        if end < start:
            raise GenBankError(f"a location that ends before it starts: {text!r}")
        return Location((Span(start - 1, end),), 1, (left_partial, right_partial))
    if "^" in text:
        left, _, right = text.partition("^")
        start, _ = _position(left)
        end, _ = _position(right)
        if end not in (start + 1, 1):
            raise GenBankError(f"a site between bases that are not adjacent: {text!r}")
        return Location((Span(start, start),), 1, (False, False))
    if "." in text:
        left, _, right = text.partition(".")
        start, _ = _position(left)
        end, _ = _position(right)
        return Location((Span(start - 1, end),), 1, (False, False))
    position, partial = _position(text)
    return Location((Span(position - 1, position),), 1, (partial, partial))


def _split(text: str) -> list[str]:
    """A comma-separated list, respecting the brackets inside it."""
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for character in text:
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        if character == "," and depth == 0:
            parts.append("".join(current))
            current = []
            continue
        current.append(character)
    parts.append("".join(current))
    return [part.strip() for part in parts if part.strip()]


def parse_location(text: str) -> Location:
    """Where a feature lies, from the location string of a GenBank record.

    Handles ``complement``, ``join`` and ``order`` to any depth, positions
    marked ``<`` or ``>``, single bases and sites written ``102^103``.
    """
    stripped = text.strip()
    if not stripped:
        raise GenBankError("a feature with no location")
    for wrapper in ("complement", "join", "order"):
        prefix = wrapper + "("
        if stripped.startswith(prefix) and stripped.endswith(")"):
            inside = stripped[len(prefix) : -1]
            if wrapper == "complement":
                within = parse_location(inside)
                return Location(
                    within.parts,
                    -within.strand,
                    (within.partial[1], within.partial[0]),
                )
            pieces = [parse_location(piece) for piece in _split(inside)]
            if not pieces:
                raise GenBankError(f"an empty {wrapper}: {text!r}")
            strands = {piece.strand for piece in pieces}
            strand = strands.pop() if len(strands) == 1 else 0
            parts = tuple(part for piece in pieces for part in piece.parts)
            return Location(
                tuple(sorted(parts)),
                strand,
                (pieces[0].partial[0], pieces[-1].partial[1]),
            )
    if stripped.startswith("one-of("):
        raise GenBankError(
            f"{text!r} names several possible positions, and a map can draw one"
        )
    return _simple(stripped)


def _unquote(value: str) -> str:
    """A qualifier value with its quotes removed and its wrapping undone.

    GenBank wraps a long value across lines at a column, so the lines are
    rejoined with a space unless the break fell inside a word, which is how
    the format writes a long translation.
    """
    if value.startswith('"') and value.endswith('"') and len(value) > 1:
        value = value[1:-1]
    return value.replace('""', '"')


#: Qualifiers whose value is a single token wrapped across lines rather than
#: prose, and so is rejoined with nothing between the pieces. A translation
#: rejoined with spaces holds characters that are not amino acids.
UNSPACED_QUALIFIERS: frozenset[str] = frozenset(
    {"translation", "peptide", "anticodon"}
)


def _qualifiers(lines: Sequence[str]) -> list[tuple[str, str]]:
    """The ``/key=value`` qualifiers of one feature, wrapping undone.

    Pairs are kept in the order the record wrote them, repeats included, so a
    feature with two ``/note`` qualifiers keeps both and everything that reads
    one takes the first.
    """
    found: list[tuple[str, str]] = []
    key: str | None = None
    pieces: list[str] = []

    def flush() -> None:
        if key is None:
            return
        separator = "" if key in UNSPACED_QUALIFIERS else " "
        joined_value = ""
        for index, piece in enumerate(pieces):
            if index and separator and not joined_value.endswith("-"):
                joined_value += separator
            joined_value += piece
        found.append((key, _unquote(joined_value.strip())))

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("/"):
            flush()
            name, separator, value = stripped[1:].partition("=")
            key = name.strip()
            pieces = [value] if separator else []
            if not separator:
                pieces = [""]
        elif key is not None:
            pieces.append(stripped)
    flush()
    return found


def _name(kind: str, qualifiers: Sequence[tuple[str, str]]) -> str:
    """The label for a feature, from the first qualifier that offers one.

    'First' means first in the record as well as first in the order of
    preference: a feature with two ``/note`` qualifiers is labelled with the
    same one that ``Feature.qualifier('note')`` returns, rather than with the
    last, which is what building a dictionary would have given.
    """
    for key in NAME_QUALIFIERS:
        for name, value in qualifiers:
            if name == key and value:
                return value.splitlines()[0].strip()
    return kind


@dataclass(frozen=True)
class Locus:
    """The LOCUS line, read field by field from the left."""

    name: str
    length: int
    molecule: str
    topology: str
    stated: bool = True

    @property
    def is_circular(self) -> bool:
        """Whether the LOCUS line said the ends are joined."""
        return self.topology == "circular"


def parse_locus(line: str) -> Locus:
    """The LOCUS line of a record.

    Read left to right: name, length, the unit, then the optional molecule
    type, topology and division. A line that stops early simply has fewer
    fields, rather than having its last field read as its topology.

    A line that names no topology is taken as linear, which is what the
    format means by leaving it out, and ``stated`` is False so that whatever
    is drawn from it can say the record did not actually say.
    """
    if not line.startswith("LOCUS"):
        raise GenBankError(f"not a LOCUS line: {line.strip()!r}")
    fields = line[len("LOCUS") :].split()
    if len(fields) < 3 or not fields[1].isdigit():
        raise GenBankError(f"a LOCUS line with no length: {line.strip()!r}")
    if fields[2] not in ("bp", "aa", "rc"):
        raise GenBankError(f"a LOCUS line with no unit: {line.strip()!r}")
    if fields[2] == "aa":
        raise GenBankError("this is a protein record, and a plasmid map is of DNA")
    rest = fields[3:]
    molecule = ""
    topology = ""
    if rest and rest[0] not in ("linear", "circular"):
        molecule = rest.pop(0)
    if rest and rest[0] in ("linear", "circular"):
        topology = rest.pop(0)
    return Locus(
        name=fields[0],
        length=int(fields[1]),
        molecule=molecule,
        topology=topology or "linear",
        stated=bool(topology),
    )


def _keyword(lines: Sequence[str], keyword: str) -> str:
    """One header field, its continuation lines folded into it."""
    collected: list[str] = []
    taking = False
    for line in lines:
        if line.startswith(keyword + " ") or line == keyword:
            taking = True
            collected.append(line[len(keyword) :].strip())
            continue
        if taking:
            if line[:1].isspace() and line.strip():
                collected.append(line.strip())
                continue
            break
    return " ".join(part for part in collected if part).strip()


def _sequence(lines: Sequence[str]) -> str:
    """The bases of an ORIGIN block, without the position numbers."""
    letters: list[str] = []
    for line in lines:
        for piece in line.split()[1:]:
            letters.append(piece)
    return "".join(letters).upper()


def split_records(text: str) -> list[list[str]]:
    """A file's records, each as its lines, split on the ``//`` terminator."""
    records: list[list[str]] = []
    current: list[str] = []
    for line in text.splitlines():
        if line.startswith("//"):
            if current:
                records.append(current)
            current = []
            continue
        current.append(line)
    if any(line.strip() for line in current):
        records.append(current)
    return records


def _features(lines: Sequence[str], locus: Locus) -> Iterator[Feature]:
    """Every feature of one record, in the order the record wrote them."""
    inside = False
    key = ""
    location: list[str] = []
    qualifier_lines: list[str] = []

    def build() -> Iterator[Feature]:
        if not key:
            return
        text = "".join(location)
        try:
            where = parse_location(text)
        except GenBankError as error:
            raise GenBankError(f"{key} at {text!r}: {error}") from None
        pairs = _qualifiers(qualifier_lines)
        try:
            parts = joined(where.parts, locus.length, locus.is_circular)
        except FeatureError as error:
            raise GenBankError(f"{key} at {text!r}: {error}") from None
        yield Feature(
            kind=key,
            name=_name(key, pairs),
            parts=parts,
            strand=where.strand,
            qualifiers=tuple(pairs),
            partial=where.partial,
        )

    for line in lines:
        if line.startswith("FEATURES"):
            inside = True
            continue
        if not inside:
            continue
        if line[:1] not in ("", " ") or line.startswith("ORIGIN"):
            break
        if not line.strip():
            continue
        if line[KEY_COLUMN : KEY_COLUMN + 1].strip() and not line.strip().startswith(
            "/"
        ):
            yield from build()
            key = line[KEY_COLUMN:LOCATION_COLUMN].strip()
            location = [line[LOCATION_COLUMN:].strip()]
            qualifier_lines = []
            continue
        if not key:
            continue
        stripped = line.strip()
        if stripped.startswith("/") or qualifier_lines:
            qualifier_lines.append(line)
        else:
            location.append(stripped)
    yield from build()


def parse_record(lines: Sequence[str]) -> Construct:
    """One GenBank record as a construct ready to draw."""
    header = next((line for line in lines if line.startswith("LOCUS")), None)
    if header is None:
        raise GenBankError("no LOCUS line, so this is not a GenBank record")
    locus = parse_locus(header)
    origin = next(
        (index for index, line in enumerate(lines) if line.startswith("ORIGIN")), None
    )
    sequence = _sequence(lines[origin + 1 :]) if origin is not None else ""
    if sequence and len(sequence) != locus.length:
        raise GenBankError(
            f"{locus.name} declares {locus.length:,} bases and its ORIGIN block "
            f"holds {len(sequence):,}"
        )
    definition = _keyword(lines, "DEFINITION").rstrip(".")
    accession = _keyword(lines, "ACCESSION").split(" ")[0]
    source = _keyword(lines, "SOURCE")
    features = tuple(_features(lines, locus))
    try:
        return Construct(
            name=locus.name,
            length=locus.length,
            topology="circular" if locus.is_circular else "linear",
            features=features,
            sequence=sequence or None,
            description=definition,
            accession=accession,
            molecule=locus.molecule,
            source=source,
        )
    except FeatureError as error:
        raise GenBankError(str(error)) from None


def parse(text: str) -> list[Construct]:
    """Every record in a GenBank file."""
    records = split_records(text)
    if not records:
        raise GenBankError("the file holds no records")
    return [parse_record(record) for record in records]


def stated_topology(lines: Sequence[str]) -> bool:
    """Whether one record's LOCUS line actually named a topology.

    A record that did not name one is drawn as linear, which is what the
    format means, but a drawing made from it should say that the record was
    silent rather than that it said linear.
    """
    header = next((line for line in lines if line.startswith("LOCUS")), None)
    if header is None:
        return False
    return parse_locus(header).stated


def looks_like(text: str) -> bool:
    """Whether a file begins the way a GenBank flat file does."""
    for line in text.splitlines():
        if not line.strip():
            continue
        return bool(re.match(r"^LOCUS\s", line))
    return False
