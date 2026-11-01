"""Getting a construct out of whatever was handed to the program.

Three ways in: a GenBank flat file, which carries the sequence, the features
and the topology; a FASTA file, which carries the sequence alone; and a
feature table or a set of command-line annotations, which carry neither
sequence nor topology. The one thing they disagree about is whether the ends
are joined, and that decides whether a feature written across the origin is
one feature or two, so it is never guessed silently: a file that states its
topology is believed, and a file that cannot state one is drawn in the shape
the command implies and told about in the summary.

A gzipped file is read as one whatever it is called, since a record downloaded
from a database usually arrives that way, and ``-`` is standard input.
"""

from __future__ import annotations

import gzip
import io
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from plasmid_map import genbank, tables
from plasmid_map.feature import (
    TOPOLOGIES,
    Construct,
    Feature,
    FeatureError,
    Span,
    joined,
)

#: The input formats, with the one that works it out from the first record.
FORMATS: tuple[str, ...] = ("genbank", "fasta")

#: What ``--topology`` takes.
TOPOLOGY_CHOICES: tuple[str, ...] = ("auto",) + TOPOLOGIES

#: The first two bytes of a gzipped file.
GZIP_MAGIC: bytes = b"\x1f\x8b"


class ReadError(ValueError):
    """An input that cannot be read."""


def read_text(path: str) -> str:
    """The whole of an input, from a file, a gzipped file or standard input."""
    if path == "-":
        data = sys.stdin.buffer.read()
    else:
        where = Path(path)
        if not where.exists():
            raise ReadError(f"there is no file at {path}")
        if where.is_dir():
            raise ReadError(f"{path} is a directory")
        data = where.read_bytes()
    if data[:2] == GZIP_MAGIC:
        try:
            data = gzip.decompress(data)
        except OSError as error:
            raise ReadError(f"{path} begins like a gzipped file but is not: {error}")
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("latin-1")


def parse_fasta(text: str) -> list[tuple[str, str, str]]:
    """Every record of a FASTA file as its name, description and sequence."""
    records: list[tuple[str, str, str]] = []
    name = ""
    description = ""
    pieces: list[str] = []
    started = False
    for line in text.splitlines():
        if line.startswith(">"):
            if started:
                records.append((name, description, "".join(pieces)))
            started = True
            heading = line[1:].strip()
            name, _space, description = heading.partition(" ")
            pieces = []
            continue
        if not started:
            if line.strip():
                raise ReadError(
                    "a FASTA file begins with a '>' line, and this one begins "
                    f"with {line.strip()[:40]!r}"
                )
            continue
        pieces.append("".join(line.split()))
    if started:
        records.append((name, description, "".join(pieces)))
    if not records:
        raise ReadError("the file holds no FASTA records")
    return records


def detect(text: str) -> str:
    """Which format an input is in, from what its first line looks like."""
    if genbank.looks_like(text):
        return "genbank"
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            return "fasta"
        break
    raise ReadError(
        "this is neither a GenBank flat file, which begins with a LOCUS "
        "line, nor a FASTA file, which begins with a '>' line; name the "
        "format with --input-format if it is one of them"
    )


@dataclass(frozen=True)
class Loaded:
    """What came out of an input, and what the input did and did not say."""

    constructs: tuple[Construct, ...]
    format: str
    topology_stated: bool
    path: str

    @property
    def first(self) -> Construct:
        """The record a map is drawn of."""
        return self.constructs[0]

    @property
    def extra(self) -> int:
        """How many further records the file held."""
        return max(len(self.constructs) - 1, 0)


def load(
    path: str, input_format: str = "auto", topology: str = "auto"
) -> Loaded:
    """Every construct in an input, in whichever format it is in."""
    if input_format not in FORMATS + ("auto",):
        raise ReadError(
            f"not an input format: {input_format!r}; "
            f"it is one of {', '.join(FORMATS)}"
        )
    if topology not in TOPOLOGY_CHOICES:
        raise ReadError(
            f"not a topology: {topology!r}; "
            f"it is one of {', '.join(TOPOLOGY_CHOICES)}"
        )
    text = read_text(path)
    if not text.strip():
        raise ReadError(f"{path} is empty")
    chosen = detect(text) if input_format == "auto" else input_format
    if chosen == "genbank":
        try:
            records = genbank.split_records(text)
            constructs = tuple(genbank.parse_record(record) for record in records)
            stated = all(genbank.stated_topology(record) for record in records)
        except genbank.GenBankError as error:
            raise ReadError(f"{path}: {error}") from None
        if not constructs:
            raise ReadError(f"{path}: the file holds no records")
    else:
        constructs = tuple(
            _from_fasta(name, description, sequence)
            for name, description, sequence in parse_fasta(text)
        )
        stated = False
    if topology != "auto":
        constructs = tuple(
            _with_topology(construct, topology) for construct in constructs
        )
        stated = True
    return Loaded(
        constructs=constructs, format=chosen, topology_stated=stated, path=path
    )


def blank(length: int, topology: str = "linear", name: str = "construct") -> Loaded:
    """A construct of a stated length with no sequence and no features.

    This is what a feature table or a set of command-line annotations is drawn
    on when there is no sequence file to draw them on: the table says where
    things are, and the length says how far round or along they go. Nothing
    here can measure GC, and asking for a GC track is refused rather than
    drawn empty.
    """
    if length <= 0:
        raise ReadError(f"a construct needs a length of at least 1, not {length}")
    if topology not in TOPOLOGIES:
        raise ReadError(
            f"not a topology: {topology!r}; it is one of {', '.join(TOPOLOGIES)}"
        )
    return Loaded(
        constructs=(Construct(name=name, length=length, topology=topology),),
        format="length",
        topology_stated=False,
        path="",
    )


def _from_fasta(name: str, description: str, sequence: str) -> Construct:
    """One FASTA record as a construct with no features and no topology."""
    if not sequence:
        raise ReadError(f"{name or 'the record'} holds no sequence")
    return Construct(
        name=name or "sequence",
        length=len(sequence),
        topology="linear",
        sequence=sequence.upper(),
        description=description,
    )


def _with_topology(construct: Construct, topology: str) -> Construct:
    """The same construct said to be circular or linear.

    Changing a circular construct to a linear one has to undo any feature that
    was read as crossing the origin, since on a linear molecule it is two
    pieces rather than one.
    """
    if construct.topology == topology:
        return construct
    circular = topology == "circular"
    features = []
    for feature in construct.features:
        parts: list[Span] = []
        for part in feature.parts:
            if part.end > construct.length:
                parts.append(Span(part.start, construct.length))
                parts.append(Span(0, part.end - construct.length))
            else:
                parts.append(part)
        # Going the other way, parts that meet at the origin are one piece on
        # a circular molecule, so the same normalisation the GenBank reader
        # does has to be done here or the override would leave a feature that
        # is two pieces on a construct where it is one.
        features.append(
            feature.with_parts(joined(tuple(sorted(parts)), construct.length, circular))
        )
    return Construct(
        name=construct.name,
        length=construct.length,
        topology=topology,
        features=tuple(features),
        sequence=construct.sequence,
        description=construct.description,
        accession=construct.accession,
        molecule=construct.molecule,
        source=construct.source,
    )


def annotate(
    construct: Construct,
    table: str = "",
    table_format: str = "auto",
    typed: Sequence[str] = (),
    seqid: str = "",
) -> tuple[Construct, int]:
    """A construct with the features of a table and the command line added.

    The added features are appended rather than replacing what the record
    already holds, so a GenBank file can be marked up with the few things
    somebody added to a vector without losing what the file already said.
    """
    added: list[Feature] = []
    if table:
        text = read_text(table)
        try:
            added.extend(
                tables.read(
                    text,
                    construct.length,
                    construct.is_circular,
                    table_format,
                    seqid,
                )
            )
        except tables.TableError as error:
            raise ReadError(f"{table}: {error}") from None
    if typed:
        try:
            added.extend(
                tables.annotations(typed, construct.length, construct.is_circular)
            )
        except tables.TableError as error:
            raise ReadError(str(error)) from None
    if not added:
        return construct, 0
    try:
        return construct.with_features(tuple(construct.features) + tuple(added)), len(
            added
        )
    except FeatureError as error:
        raise ReadError(str(error)) from None


def open_output(path: str) -> io.TextIOBase:
    """A file to write to, or standard output for ``-``."""
    if path == "-":
        return sys.stdout
    try:
        return open(path, "w", encoding="utf-8", newline="\n")
    except OSError as error:
        raise ReadError(f"cannot write to {path}: {error}") from None
