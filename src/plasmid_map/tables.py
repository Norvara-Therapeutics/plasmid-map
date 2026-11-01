"""Features written as a table, or typed on the command line.

Not every construct has a GenBank file. A cassette being designed, a synthesis
order, a region someone measured: these are a name and two numbers, and this
module turns them into features.

Two table formats are read. A delimited table - tab, comma or semicolon, with
a header row - counts from 1 and includes both ends, because that is how a
person writing a table by hand writes one and how every sequence viewer shows
a position. GFF3 counts the same way and is read as itself. Which convention a
table was read in is never guessed from the numbers; it follows from the
format, and the format follows from the header.

A table may give positions as ``start`` and ``end`` columns or as a single
``location`` column written the way GenBank writes one, which is how a spliced
or origin-crossing feature is put in a table at all.

The two say different things about the strand, and deliberately. A location
written the GenBank way carries GenBank's meaning, where a location without
``complement()`` is on the plus strand; a ``start`` and an ``end`` column say
nothing about strand at all unless a strand column does, and a feature with no
strand is drawn without a point on it. A strand column, where there is one,
settles it either way.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from typing import Iterator, Sequence

from plasmid_map.feature import Feature, FeatureError, Span, joined, strand_of
from plasmid_map.genbank import GenBankError, parse_location

#: The column holding a feature's label, in the order the names are tried.
NAME_COLUMNS: tuple[str, ...] = ("name", "label", "feature_name", "gene", "id")

#: The column holding a feature's kind.
KIND_COLUMNS: tuple[str, ...] = ("kind", "type", "feature", "feature_type")

#: The column holding a strand.
STRAND_COLUMNS: tuple[str, ...] = ("strand", "sense", "direction")

#: The column holding a whole GenBank-style location.
LOCATION_COLUMNS: tuple[str, ...] = ("location", "position", "span")

#: The columns holding the two ends of a simple feature.
START_COLUMNS: tuple[str, ...] = ("start", "from", "begin", "left")
END_COLUMNS: tuple[str, ...] = ("end", "to", "finish", "right")

#: What a feature with no kind of its own is called.
DEFAULT_KIND: str = "misc_feature"

#: The delimiters tried when reading a table, in order.
DELIMITERS: str = "\t,;"

#: The formats this module reads.
FORMATS: tuple[str, ...] = ("table", "gff3")


class TableError(ValueError):
    """A feature table that cannot be read."""


def _pick(header: Sequence[str], names: Sequence[str]) -> int | None:
    """Which column of a header is the first of these names, if any."""
    lowered = [item.strip().lower().lstrip("#") for item in header]
    for name in names:
        if name in lowered:
            return lowered.index(name)
    return None


def _delimiter(text: str) -> str:
    """Which delimiter a table uses, from its first line."""
    first = next((line for line in text.splitlines() if line.strip()), "")
    counts = {character: first.count(character) for character in DELIMITERS}
    best = max(counts, key=lambda character: counts[character])
    if not counts[best]:
        raise TableError(
            "the first line holds no tab, comma or semicolon, so this is not "
            "a table with a header row"
        )
    return best


def _value(row: Sequence[str], index: int | None) -> str:
    if index is None or index >= len(row):
        return ""
    return row[index].strip()


def _location(text: str, length: int, circular: bool) -> tuple[tuple[Span, ...], int]:
    """A GenBank-style location from a table cell."""
    try:
        where = parse_location(text)
    except GenBankError as error:
        raise TableError(f"{text!r}: {error}") from None
    return joined(where.parts, length, circular), where.strand


def _simple(start: str, end: str) -> Span:
    """A span from two columns counting from 1 with both ends included."""
    try:
        first = int(start)
        last = int(end)
    except ValueError:
        raise TableError(
            f"a start and an end have to be whole numbers, not "
            f"{start!r} and {end!r}"
        ) from None
    if first < 1:
        raise TableError(f"a table counts from 1, so {first} is not a position")
    if last < first:
        raise TableError(f"a feature that ends before it starts: {first}..{last}")
    return Span(first - 1, last)


def read_table(text: str, length: int, circular: bool) -> list[Feature]:
    """Features from a delimited table with a header row.

    Required: a name column and either a location column or a start and an end
    column. Everything else - kind, strand, colour, any further column - is
    optional, and a column this module has no use for is carried through as a
    qualifier rather than dropped.
    """
    delimiter = _delimiter(text)
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    rows = [row for row in rows if row and any(cell.strip() for cell in row)]
    if not rows:
        raise TableError("the table holds no rows")
    header = [cell.strip().lstrip("#") for cell in rows[0]]
    name_at = _pick(header, NAME_COLUMNS)
    location_at = _pick(header, LOCATION_COLUMNS)
    start_at = _pick(header, START_COLUMNS)
    end_at = _pick(header, END_COLUMNS)
    if location_at is None and (start_at is None or end_at is None):
        raise TableError(
            "a table needs a location column, or a start and an end column; "
            f"this one has {', '.join(header) or 'no header'}"
        )
    kind_at = _pick(header, KIND_COLUMNS)
    strand_at = _pick(header, STRAND_COLUMNS)
    spoken = {name_at, location_at, start_at, end_at, kind_at, strand_at}
    features: list[Feature] = []
    for number, row in enumerate(rows[1:], start=2):
        try:
            if location_at is not None and _value(row, location_at):
                parts, strand = _location(
                    _value(row, location_at), length, circular
                )
            else:
                parts = (_simple(_value(row, start_at), _value(row, end_at)),)
                strand = 0
            written = _value(row, strand_at)
            if written:
                strand = strand_of(written)
            extra = tuple(
                (header[index], cell.strip())
                for index, cell in enumerate(row)
                if index < len(header) and index not in spoken and cell.strip()
            )
            features.append(
                Feature(
                    kind=_value(row, kind_at) or DEFAULT_KIND,
                    name=_value(row, name_at),
                    parts=parts,
                    strand=strand,
                    qualifiers=extra,
                )
            )
        except (TableError, FeatureError) as error:
            raise TableError(f"line {number}: {error}") from None
    return features


def _attributes(text: str) -> list[tuple[str, str]]:
    """The ninth column of a GFF3 line, as key and value pairs."""
    from urllib.parse import unquote

    pairs: list[tuple[str, str]] = []
    for piece in text.split(";"):
        if not piece.strip():
            continue
        key, _separator, value = piece.strip().partition("=")
        pairs.append((unquote(key.strip()), unquote(value.strip())))
    return pairs


#: Which GFF3 attribute gives a feature its label, most preferred first.
GFF_NAMES: tuple[str, ...] = ("Name", "gene", "gene_name", "label", "ID", "locus_tag")


@dataclass
class _Group:
    """The rows of a GFF3 file that are parts of one feature."""

    kind: str
    name: str
    strand: int
    qualifiers: list[tuple[str, str]]
    spans: list[Span]


#: How many tab-separated columns a GFF3 line has. Exactly nine: a line with
#: eight has no attributes, which is where the name and the grouping come
#: from, and a line with ten has a tab inside a field.
GFF_COLUMNS: int = 9


def seqids(text: str) -> list[str]:
    """Every distinct sequence name a GFF3 file names, in the order it does."""
    found: list[str] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name = line.split("\t")[0]
        if name not in found:
            found.append(name)
    return found


def read_gff3(
    text: str, length: int, circular: bool, seqid: str = ""
) -> list[Feature]:
    """Features from a GFF3 file, which counts from 1 with both ends in.

    Rows sharing an ``ID`` attribute are one feature with several parts, which
    is how GFF3 writes a spliced coding sequence; building one feature per row
    instead would lose the join and count the feature twice.

    Which sequence a file describes cannot be checked from here - a GFF3 file
    and a GenBank record may legitimately name the same molecule differently -
    so a file naming more than one sequence is refused rather than having its
    rows silently added together on one construct. Name the one to use and it
    is used.
    """
    present = seqids(text)
    if seqid and seqid not in present:
        raise TableError(
            f"the table names no sequence {seqid!r}; it names "
            f"{', '.join(repr(name) for name in present) or 'none'}"
        )
    if not seqid and len(present) > 1:
        raise TableError(
            f"the table describes {len(present)} sequences "
            f"({', '.join(repr(name) for name in present)}), and a map is of "
            "one; choose it with --features-seqid"
        )
    grouped: dict[str, _Group] = {}
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != GFF_COLUMNS:
            raise TableError(
                f"line {number}: a GFF3 line has {GFF_COLUMNS} tab-separated "
                f"columns, and this one has {len(fields)}"
            )
        if seqid and fields[0] != seqid:
            continue
        try:
            span = _simple(fields[3], fields[4])
        except TableError as error:
            raise TableError(f"line {number}: {error}") from None
        pairs = _attributes(fields[8])
        table = dict(pairs)
        name = next((table[key] for key in GFF_NAMES if table.get(key)), "")
        kind = fields[2] or DEFAULT_KIND
        strand = strand_of(fields[6])
        identity = table.get("ID", "")
        key = (
            f"{fields[0]}\t{kind}\t{strand}\t{identity}"
            if identity
            else f"line {number}"
        )
        group = grouped.get(key)
        if group is None:
            grouped[key] = _Group(
                kind=kind, name=name, strand=strand, qualifiers=pairs, spans=[span]
            )
        else:
            group.spans.append(span)
            if not group.name:
                group.name = name
    return [
        Feature(
            kind=group.kind,
            name=group.name,
            parts=joined(tuple(group.spans), length, circular),
            strand=group.strand,
            qualifiers=tuple(group.qualifiers),
        )
        for group in grouped.values()
    ]


def looks_like_gff3(text: str) -> bool:
    """Whether a file is GFF3 rather than a table with a header row.

    The pragma settles it. Without one, a nine-column line is GFF3 only if it
    reads like a GFF3 record: two whole numbers where the coordinates go and
    something that is a strand where the strand goes. A header row naming nine
    columns passes neither test.
    """
    for line in text.splitlines():
        if line.startswith("##gff-version"):
            return True
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != GFF_COLUMNS:
            return False
        if not (fields[3].strip().isdigit() and fields[4].strip().isdigit()):
            return False
        return fields[6].strip() in ("+", "-", ".", "?")
    return False


def read(
    text: str,
    length: int,
    circular: bool,
    kind: str = "auto",
    seqid: str = "",
) -> list[Feature]:
    """Features from a table, in whichever of the two formats it is in."""
    if kind not in FORMATS + ("auto",):
        raise TableError(
            f"not a feature table format: {kind!r}; "
            f"it is one of {', '.join(FORMATS)}"
        )
    chosen = kind
    if chosen == "auto":
        chosen = "gff3" if looks_like_gff3(text) else "table"
    if chosen == "gff3":
        return read_gff3(text, length, circular, seqid)
    if seqid:
        raise TableError(
            "--features-seqid names a sequence in a GFF3 file, and this table "
            "is not one"
        )
    return read_table(text, length, circular)


@dataclass(frozen=True)
class Annotation:
    """One feature typed on the command line."""

    name: str
    start: int
    end: int
    strand: int
    kind: str


def parse_annotation(text: str) -> Annotation:
    """One ``--annotate`` argument.

    Written ``NAME:START..END``, with an optional strand and kind after it:
    ``rep:1..1200:+:CDS``. A single position is written as one number. The
    name is whatever comes before the positions, so a name may hold a colon
    even though the fields are separated by one.
    """
    pieces = text.split(":")
    where = None
    for index, piece in enumerate(pieces):
        stripped = piece.strip()
        if stripped and all(
            part.strip().isdigit() for part in stripped.split("..") if part.strip()
        ):
            if ".." in stripped or stripped.isdigit():
                where = index
                break
    if where is None or where == 0:
        raise TableError(
            f"{text!r} is not an annotation; write it as NAME:START..END, "
            "with an optional strand and kind after it"
        )
    name = ":".join(pieces[:where]).strip()
    positions = pieces[where].strip()
    rest = [piece.strip() for piece in pieces[where + 1 :]]
    if ".." in positions:
        first, _separator, last = positions.partition("..")
    else:
        first = last = positions
    span = _simple(first, last)
    strand = strand_of(rest[0]) if rest and rest[0] else 0
    kind = rest[1] if len(rest) > 1 and rest[1] else DEFAULT_KIND
    return Annotation(
        name=name, start=span.start, end=span.end, strand=strand, kind=kind
    )


def annotations(texts: Sequence[str], length: int, circular: bool) -> list[Feature]:
    """Every ``--annotate`` argument as a feature."""
    built: list[Feature] = []
    for text in texts:
        found = parse_annotation(text)
        span = Span(found.start, found.end)
        if span.end > length and not circular:
            raise TableError(
                f"{found.name} ends at {span.end:,}, past the {length:,} bases "
                "of a linear construct"
            )
        built.append(
            Feature(
                kind=found.kind,
                name=found.name,
                parts=(span,),
                strand=found.strand,
            )
        )
    return built


def describe() -> Iterator[tuple[str, str]]:
    """What each readable table format is, for the help and the README."""
    yield (
        "table",
        "tab, comma or semicolon separated with a header row; counts from 1 "
        "with both ends included",
    )
    yield (
        "gff3",
        "GFF3; counts from 1 with both ends included, and takes its label "
        "from the Name attribute",
    )
