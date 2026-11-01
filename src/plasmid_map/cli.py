"""The command line.

Four commands. ``circular`` and ``linear`` draw the two shapes a construct can
be shown in; ``features`` prints what would be drawn without drawing it, which
is how a filter is checked before a figure is made; ``types`` lists the feature
kinds the program knows and the colour each is given.

The two drawing commands take the same options, because the question 'which
features, labelled how, in what colours' has the same answer whichever shape
the construct is shown in. An option that would mean nothing in one of the two
- ``--start``, which turns a circle, and has nothing to turn on a line - is
not offered there at all, so asking for it is an error rather than an option
that is quietly ignored.

Everything that can be refused is refused before anything is read, and the
outputs are checked against each other and against the input by resolved path,
so a run cannot write two of its outputs over each other or overwrite the file
it is drawing from.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import sys
from dataclasses import replace
from pathlib import Path
from typing import IO, Any, Sequence

from plasmid_map import circular, linear, reader, report
from plasmid_map import rows as tabular
from plasmid_map import text as prose
from plasmid_map.composition import CompositionError, measure
from plasmid_map.drawing import LABEL_MODES, Drawing, DrawingError, Style, select
from plasmid_map.feature import FeatureError
from plasmid_map.labels import LABEL_WIDTH
from plasmid_map.palette import PALETTES, ColourError, Palette
from plasmid_map.reader import ReadError
from plasmid_map.tables import FORMATS as TABLE_FORMATS

#: What may be written to standard output.
STDOUT_FORMATS: tuple[str, ...] = ("table", "tsv", "json", "svg", "none")

#: What ``--input-format`` takes, with the default that works it out.
INPUT_FORMATS: tuple[str, ...] = ("auto",) + reader.FORMATS

#: What ``--features-format`` takes.
FEATURE_FORMATS: tuple[str, ...] = ("auto",) + TABLE_FORMATS

#: How many rows of a table standard output shows before pointing elsewhere.
PREVIEW_ROWS: int = 25


class UsageError(ValueError):
    """A combination of options that cannot be carried out."""


def _positive(value: str) -> int:
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a whole number: {value!r}") from None
    if number < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1, not {value!r}")
    return number


def _non_negative(value: str) -> int:
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a whole number: {value!r}") from None
    if number < 0:
        raise argparse.ArgumentTypeError(f"cannot be negative: {value!r}")
    return number


def _measure(value: str) -> float:
    """A length on the canvas, in the units the SVG is drawn in."""
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {value!r}") from None
    if number != number or number in (float("inf"), float("-inf")):
        raise argparse.ArgumentTypeError(f"not a length: {value!r}")
    if number <= 0:
        raise argparse.ArgumentTypeError(f"must be more than zero, not {value!r}")
    return number


def _kinds(value: str) -> list[str]:
    """Feature kinds, comma-separated, as a GenBank file writes them."""
    names = [part.strip() for part in value.split(",") if part.strip()]
    if not names:
        raise argparse.ArgumentTypeError("name at least one kind")
    return names


def _colour(value: str) -> tuple[str, str]:
    """A ``KIND=COLOUR`` pair for overriding one kind's colour."""
    kind, separator, colour = value.partition("=")
    if not separator or not kind.strip() or not colour.strip():
        raise argparse.ArgumentTypeError(
            f"write a colour override as KIND=COLOUR, not {value!r}"
        )
    from plasmid_map.palette import check_colour

    try:
        return kind.strip(), check_colour(colour)
    except ColourError as error:
        raise argparse.ArgumentTypeError(str(error)) from None


def _version() -> str:
    from plasmid_map import __version__

    return __version__


def _add_input(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "input",
        metavar="INPUT",
        nargs="?",
        help=(
            "a GenBank or FASTA file; '-' for standard input, gzip read as "
            "read. Leave it out and give --length to draw a feature table or "
            "command-line annotations on their own"
        ),
    )
    parser.add_argument(
        "--length",
        metavar="BASES",
        type=_positive,
        help=(
            "how long the construct is, for drawing a feature table or "
            "annotations with no sequence file; not for use with an input"
        ),
    )
    parser.add_argument(
        "--input-format",
        choices=INPUT_FORMATS,
        default="auto",
        help=(
            "the format of the input; by default it is worked out from its "
            "first line (default: auto)"
        ),
    )
    parser.add_argument(
        "--topology",
        choices=reader.TOPOLOGY_CHOICES,
        default="auto",
        help=(
            "whether the ends are joined; by default a GenBank record is "
            "believed and anything else is taken as linear (default: auto)"
        ),
    )
    parser.add_argument(
        "--features",
        metavar="PATH",
        help=(
            "a feature table to add to whatever the input already carries, "
            "as a delimited table with a header row or as GFF3"
        ),
    )
    parser.add_argument(
        "--features-format",
        choices=FEATURE_FORMATS,
        default="auto",
        help="the format of --features (default: auto)",
    )
    parser.add_argument(
        "--features-seqid",
        metavar="NAME",
        default="",
        help=(
            "which sequence of a GFF3 feature table to use, named as its "
            "first column names it; needed when the file describes more than "
            "one, since a map is of one"
        ),
    )
    parser.add_argument(
        "--record",
        metavar="NAME-OR-NUMBER",
        default="",
        help=(
            "which record of a file holding several to draw, by accession, by "
            "name, or by its position counting from 1 (default: the first)"
        ),
    )
    parser.add_argument(
        "--annotate",
        metavar="NAME:START..END",
        action="append",
        default=[],
        help=(
            "one more feature, counting from 1 with both ends included, with "
            "an optional strand and kind after it, as in "
            "'EGFP:1200..1917:+:CDS'; repeatable"
        ),
    )


def _add_selection(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--kind",
        metavar="TYPE[,TYPE]",
        type=_kinds,
        default=[],
        help=(
            "draw only features of these kinds, as a GenBank file writes them "
            "in its feature table, for instance CDS or rep_origin"
        ),
    )
    parser.add_argument(
        "--hide",
        metavar="TYPE[,TYPE]",
        type=_kinds,
        default=[],
        help="draw everything except features of these kinds",
    )
    parser.add_argument(
        "--minimum-length",
        metavar="BASES",
        type=_non_negative,
        default=0,
        help="leave out features shorter than this (default: 0, leave none out)",
    )
    parser.add_argument(
        "--qualifiers",
        metavar="KEY[,KEY]",
        type=_kinds,
        default=[],
        help=(
            "carry these qualifiers of each feature across as extra columns "
            "of the tables, for instance product or note; a feature without "
            "one gets an empty field"
        ),
    )


def _add_style(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--width",
        metavar="PIXELS",
        type=_measure,
        default=760.0,
        help="how wide the drawing is, in the units the SVG uses (default: 760)",
    )
    parser.add_argument(
        "--track",
        metavar="PIXELS",
        type=_measure,
        default=15.0,
        help=(
            "how thick a feature band is; it is thinned automatically when a "
            "construct needs more lanes than will fit (default: 15)"
        ),
    )
    parser.add_argument(
        "--labels",
        choices=LABEL_MODES,
        default="outside",
        help="where the labels go, or 'none' for a map without them",
    )
    parser.add_argument(
        "--label-width",
        metavar="CHARACTERS",
        type=_positive,
        default=LABEL_WIDTH,
        help=(
            f"how many characters a label keeps before it is cut; the whole "
            f"text stays in the drawing as a tooltip (default: {LABEL_WIDTH})"
        ),
    )
    parser.add_argument(
        "--palette",
        choices=PALETTES,
        default="colour",
        help="'colour' or 'mono', which is the same map in greys (default: colour)",
    )
    parser.add_argument(
        "--colour",
        metavar="KIND=COLOUR",
        type=_colour,
        action="append",
        default=[],
        help="give one kind of feature a colour of your own; repeatable",
    )
    parser.add_argument(
        "--title",
        metavar="TEXT",
        default="",
        help="what to call the construct, instead of the name in the record",
    )
    parser.add_argument(
        "--no-ruler",
        dest="ruler",
        action="store_false",
        help="leave off the ruler and its tick labels",
    )
    parser.add_argument(
        "--legend",
        action="store_true",
        help="add a key naming the colour of each group of feature kinds",
    )
    parser.add_argument(
        "--gc",
        action="store_true",
        help="add a GC content track, which needs the sequence",
    )
    parser.add_argument(
        "--skew",
        action="store_true",
        help="add a GC skew track, (G - C) / (G + C), which needs the sequence",
    )
    parser.add_argument(
        "--gc-window",
        metavar="BASES",
        type=_positive,
        help=(
            "how wide a GC window is; by default a fortieth of the construct, "
            "never finer than 50 bases"
        ),
    )
    parser.add_argument(
        "--gc-step",
        metavar="BASES",
        type=_positive,
        help="how far apart the GC windows are placed (default: a fifth of the window)",
    )


def _add_output(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-f",
        "--format",
        choices=STDOUT_FORMATS,
        default="table",
        help=(
            "what goes to standard output: a summary as a table, the features "
            "as tsv or json, the drawing itself as svg, or nothing "
            "(default: table)"
        ),
    )
    parser.add_argument(
        "--svg",
        metavar="PATH",
        help="the drawing, as a standalone SVG file; '-' for standard output",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="the drawn features, as a tab-separated table",
    )
    parser.add_argument(
        "--report",
        metavar="PATH",
        help="a standalone HTML view: the map, what it says, and its tables",
    )


def build_parser() -> argparse.ArgumentParser:
    """The whole command line, as a parser."""
    parser = argparse.ArgumentParser(
        prog="plasmid-map",
        description=(
            "Draw publication-ready circular plasmid maps and linear feature "
            "tracks as standalone SVG, from GenBank files, feature tables or "
            "annotations given on the command line."
        ),
        epilog=(
            "Positions given to and printed by this program count from 1 and "
            "include both ends, as GenBank does. A feature crossing the origin "
            "of a circular construct is one feature, and its length is the sum "
            "of its parts."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"plasmid-map {_version()}"
    )
    commands = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    for name, description in (
        (
            "circular",
            "Draw a construct as a circle, with the features on concentric "
            "rings and their labels round the outside.",
        ),
        (
            "linear",
            "Draw a construct as a horizontal track, with the features in "
            "stacked rows and their labels above them.",
        ),
    ):
        drawer = commands.add_parser(
            name,
            help=f"draw a {name} map",
            description=description,
        )
        _add_input(drawer)
        _add_selection(drawer)
        _add_style(drawer)
        _add_output(drawer)
        if name == "circular":
            drawer.add_argument(
                "--start",
                metavar="BASE",
                type=_positive,
                default=1,
                help=(
                    "which base to put at twelve o'clock, counting from 1, "
                    "which is how a map is turned so that a feature of "
                    "interest sits at the top (default: 1)"
                ),
            )

    listing = commands.add_parser(
        "features",
        help="list the features that would be drawn, without drawing them",
        description=(
            "Print the features a map would show, after the filters have been "
            "applied and the lanes assigned, which is how a selection is "
            "checked before a figure is made."
        ),
    )
    _add_input(listing)
    _add_selection(listing)
    listing.add_argument(
        "-f",
        "--format",
        choices=("table", "tsv", "json"),
        default="table",
        help="what goes to standard output (default: table)",
    )
    listing.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="the features, as a tab-separated table",
    )
    listing.add_argument(
        "--shape",
        choices=("circular", "linear"),
        default="circular",
        help=(
            "which map's lanes to work out, since the two assign them the "
            "same way but a linear one cuts at the origin (default: circular)"
        ),
    )

    kinds = commands.add_parser(
        "types",
        help="list the feature kinds this program knows and their colours",
        description=(
            "Every feature kind with a colour of its own, the group it belongs "
            "to and the colour it is drawn in. Anything not listed is drawn in "
            "the default colour for its group."
        ),
    )
    kinds.add_argument(
        "--palette",
        choices=PALETTES,
        default="colour",
        help="which palette to list (default: colour)",
    )
    kinds.add_argument(
        "--colour",
        metavar="KIND=COLOUR",
        type=_colour,
        action="append",
        default=[],
        help="show the listing with one kind's colour overridden; repeatable",
    )
    kinds.add_argument(
        "-f",
        "--format",
        choices=("table", "tsv", "json"),
        default="table",
        help="what goes to standard output (default: table)",
    )

    return parser


def _style(arguments: argparse.Namespace) -> Style:
    """Everything about how a map should look, from the options given."""
    palette = Palette(arguments.palette, dict(arguments.colour))
    return Style(
        palette=palette,
        width=arguments.width,
        track=arguments.track,
        labels=arguments.labels,
        label_width=arguments.label_width,
        ruler=arguments.ruler,
        gc=arguments.gc,
        skew=arguments.skew,
        legend=arguments.legend,
        title=arguments.title,
        start=getattr(arguments, "start", 1) - 1,
    )


def _check_stdout(arguments: argparse.Namespace) -> None:
    """Refuse two things being written into one stream.

    An SVG document followed by a JSON document in the same stream is neither,
    and nothing downstream can take it apart again, so writing both is refused
    here rather than produced and left to be noticed.
    """
    writers: list[str] = []
    for option in ("svg", "report", "output"):
        if getattr(arguments, option, None) == "-":
            writers.append(f"--{option}" if option != "output" else "--output")
    if getattr(arguments, "format", "none") != "none":
        writers.append(f"--format {arguments.format}")
    if len(writers) > 1:
        raise UsageError(
            f"{' and '.join(writers)} all write to standard output, and one "
            "stream cannot hold both; send one to a file, or silence the "
            "summary with -f none"
        )


#: Options that only mean something alongside another option. Asking for one
#: without the other is refused rather than quietly ignored, which is the rule
#: the rest of this program follows.
DEPENDENT: tuple[tuple[str, tuple[str, ...], str], ...] = (
    ("gc_window", ("gc", "skew"), "--gc-window sizes a track; add --gc or --skew"),
    ("gc_step", ("gc", "skew"), "--gc-step spaces a track; add --gc or --skew"),
    (
        "features_seqid",
        ("features",),
        "--features-seqid chooses a sequence of a feature table; add --features",
    ),
)


def _check_options(arguments: argparse.Namespace) -> None:
    """Refuse options that would have no effect as things stand."""
    for name, needs, message in DEPENDENT:
        if getattr(arguments, name, None) and not any(
            getattr(arguments, other, None) for other in needs
        ):
            raise UsageError(message)
    if getattr(arguments, "features_format", "auto") != "auto" and not getattr(
        arguments, "features", None
    ):
        raise UsageError("--features-format names the format of --features; add it")
    if arguments.input is None:
        if getattr(arguments, "input_format", "auto") != "auto":
            raise UsageError(
                "--input-format names the format of an input, and there is no "
                "input to name; drop it, or give a file to draw"
            )
        if getattr(arguments, "record", ""):
            raise UsageError(
                "--record chooses one record of a file, and there is no file; "
                "drop it, or give one to draw"
            )
    if (
        getattr(arguments, "labels", "outside") == "none"
        and getattr(arguments, "label_width", LABEL_WIDTH) != LABEL_WIDTH
    ):
        raise UsageError(
            "--label-width cuts labels, and --labels none draws none; drop one"
        )


def _check_paths(arguments: argparse.Namespace) -> None:
    """Refuse outputs that would overwrite the input or one another."""
    _check_stdout(arguments)
    _check_options(arguments)
    named: list[tuple[str, str]] = []
    for option in ("svg", "output", "report"):
        value = getattr(arguments, option, None)
        if value and value != "-":
            named.append((f"--{option}", value))
    resolved: dict[Path, str] = {}
    for option, value in named:
        where = Path(value).resolve()
        if where in resolved:
            raise UsageError(
                f"{option} and {resolved[where]} both write to {value}; "
                "give them different paths"
            )
        resolved[where] = option
    if arguments.input and arguments.input != "-":
        source = Path(arguments.input).resolve()
        for option, value in named:
            if Path(value).resolve() == source:
                raise UsageError(
                    f"{option} would write over {arguments.input}, which is "
                    "the file being read"
                )
    if getattr(arguments, "features", None):
        table = Path(arguments.features).resolve()
        for option, value in named:
            if Path(value).resolve() == table:
                raise UsageError(
                    f"{option} would write over {arguments.features}, which is "
                    "the feature table being read"
                )


def _source(arguments: argparse.Namespace) -> str:
    """What to call the input in a message, when there may not be one."""
    return arguments.input or "the command line"


def _load(arguments: argparse.Namespace):
    """The construct to draw, with whatever was added to it."""
    if arguments.input is None:
        if not arguments.length:
            raise UsageError(
                "give a GenBank or FASTA file to draw, or --length and a "
                "feature table or --annotate to draw without one"
            )
        if not (arguments.features or arguments.annotate):
            raise UsageError(
                "--length says how long the construct is, but nothing says "
                "what is on it; add --features or --annotate"
            )
        loaded = reader.blank(
            arguments.length,
            "linear" if arguments.topology == "auto" else arguments.topology,
            getattr(arguments, "title", "") or "construct",
        )
    else:
        if arguments.length:
            raise UsageError(
                f"--length is for drawing without an input, and "
                f"{arguments.input} states its own length"
            )
        loaded = reader.load(
            arguments.input, arguments.input_format, arguments.topology
        )
    construct = _chosen(loaded, arguments.record)
    construct, added = reader.annotate(
        construct,
        table=getattr(arguments, "features", "") or "",
        table_format=getattr(arguments, "features_format", "auto"),
        typed=arguments.annotate,
        seqid=getattr(arguments, "features_seqid", "") or "",
    )
    return loaded, construct, added


def _chosen(loaded, asked: str):
    """Which record of the file to draw.

    By accession, by name, or by position counting from 1. Without ``--record``
    the first is drawn, and the summary says how many others the file held.
    """
    if not asked:
        return loaded.first
    if asked.isdigit():
        position = int(asked)
        if not 1 <= position <= len(loaded.constructs):
            raise UsageError(
                f"--record {asked} asks for record {position} of a file "
                f"holding {len(loaded.constructs):,}"
            )
        return loaded.constructs[position - 1]
    for construct in loaded.constructs:
        if asked in (construct.name, construct.accession):
            return construct
    named = ", ".join(
        sorted({item.accession or item.name for item in loaded.constructs})
    )
    raise UsageError(f"the file holds no record called {asked!r}; it holds {named}")


def _profile(construct, arguments: argparse.Namespace, shape: str):
    """The GC profile, when one was asked for and can be measured.

    Measured the way the map is drawn, not the way the molecule is: a circular
    construct shown as a line has been cut, so its windows stop at the cut
    too. Measuring round the origin and drawing on a line would put three
    windows on the same rightmost point, which is a spike that is not in the
    sequence.
    """
    if not (arguments.gc or arguments.skew):
        return None
    if construct.sequence is None:
        raise UsageError(
            "a GC track needs the sequence, and the input carries none; "
            "draw from a GenBank or FASTA file, or leave off --gc and --skew"
        )
    return measure(
        construct.sequence,
        construct.is_circular and shape == "circular",
        arguments.gc_window,
        arguments.gc_step,
    )


def _write(path: str, body: str, stdout: IO[str]) -> None:
    """One output, to a file or to standard output."""
    if path == "-":
        stdout.write(body)
        return
    with contextlib.closing(reader.open_output(path)) as handle:
        handle.write(body)


def _json_payload(
    drawing: Drawing,
    source: str,
    source_format: str,
    records: int = 1,
    stated: bool = True,
    qualifiers: Sequence[str] = (),
) -> dict[str, Any]:
    """Everything the summary says, as a structure rather than as a table."""
    columns = list(tabular.feature_columns(qualifiers))
    body = tabular.feature_rows(
        drawing.layout.placements, drawing.construct.length, qualifiers
    )
    return {
        "source": source,
        "format": source_format,
        "shape": drawing.shape,
        "construct": {
            "name": drawing.construct.name,
            "accession": drawing.construct.accession,
            "length": drawing.construct.length,
            "topology": drawing.construct.topology,
            "description": drawing.construct.description,
        },
        "canvas": {"width": drawing.width, "height": drawing.height},
        "summary": dict(tabular.summary_rows(drawing)),
        "notes": prose.notes(drawing, records, stated),
        "features": [dict(zip(columns, row)) for row in body],
    }


def _preview(columns: Sequence[str], body: Sequence[Sequence[object]]) -> str:
    """A table for a terminal, cut off before it fills the screen."""
    shown = list(body[:PREVIEW_ROWS])
    rendered = prose.table(list(columns), shown, numeric=(2, 3, 4, 6, 9, 10), indent="  ")
    if len(body) > PREVIEW_ROWS:
        rendered += (
            f"\n  \u2026 and {len(body) - PREVIEW_ROWS:,} more; "
            "write them all with --output"
        )
    return rendered


def run_draw(
    arguments: argparse.Namespace, stdout: IO[str], shape: str
) -> int:
    """Draw a map and write whatever was asked for."""
    _check_paths(arguments)
    loaded, construct, added = _load(arguments)
    profile = _profile(construct, arguments, shape)
    selection = select(
        construct,
        kinds=arguments.kind,
        hide=arguments.hide,
        minimum_length=arguments.minimum_length,
    )
    if not selection.kept:
        raise UsageError(
            f"{_source(arguments)} has nothing left to draw once the filters have "
            f"been applied; {len(construct.features):,} features were read"
        )
    style = _style(arguments)
    draw = circular.draw if shape == "circular" else linear.draw
    drawing = replace(draw(construct, selection, style, profile), added=added)
    columns = tabular.feature_columns(arguments.qualifiers)

    if arguments.svg:
        _write(arguments.svg, drawing.svg, stdout)
    if arguments.report:
        _write(
            arguments.report,
            report.build(
                drawing,
                _source(arguments),
                loaded.format,
                len(loaded.constructs),
                loaded.topology_stated,
            ),
            stdout,
        )
    body = tabular.feature_rows(
        drawing.layout.placements, construct.length, arguments.qualifiers
    )
    if arguments.output:
        _write(arguments.output, tabular.render(columns, body), stdout)

    if arguments.format == "svg":
        stdout.write(drawing.svg)
    elif arguments.format == "tsv":
        stdout.write(tabular.render(columns, body))
    elif arguments.format == "json":
        stdout.write(
            json.dumps(
                _json_payload(
                    drawing,
                    _source(arguments),
                    loaded.format,
                    len(loaded.constructs),
                    loaded.topology_stated,
                    arguments.qualifiers,
                ),
                indent=2,
            )
            + "\n"
        )
    elif arguments.format == "table":
        stdout.write(
            prose.summary(
                drawing,
                _source(arguments),
                loaded.format,
                len(loaded.constructs),
                loaded.topology_stated,
            )
        )
    return 0


def run_features(arguments: argparse.Namespace, stdout: IO[str]) -> int:
    """List what would be drawn, without drawing it."""
    _check_paths(arguments)
    loaded, construct, _added = _load(arguments)
    selection = select(
        construct,
        kinds=arguments.kind,
        hide=arguments.hide,
        minimum_length=arguments.minimum_length,
    )
    from plasmid_map.feature import minimum_width
    from plasmid_map.lanes import assign

    layout = assign(
        selection.kept,
        construct.length,
        construct.is_circular and arguments.shape == "circular",
        minimum_width(construct.length),
    )
    columns = list(tabular.feature_columns(arguments.qualifiers))
    body = tabular.feature_rows(
        layout.placements, construct.length, arguments.qualifiers
    )
    if arguments.output:
        _write(arguments.output, tabular.render(columns, body), stdout)
    if arguments.format == "tsv":
        stdout.write(tabular.render(columns, body))
        return 0
    if arguments.format == "json":
        stdout.write(
            json.dumps(
                {
                    "source": _source(arguments),
                    "construct": construct.name,
                    "length": construct.length,
                    "topology": construct.topology,
                    "lanes": layout.lanes,
                    "features": [dict(zip(columns, row)) for row in body],
                },
                indent=2,
            )
            + "\n"
        )
        return 0
    stdout.write(
        f"{_source(arguments)} ({loaded.format}, {construct.topology}, "
        f"{construct.length:,} bp)\n"
    )
    stdout.write(
        prose.definitions(
            [
                ("features in the record", f"{len(construct.features):,}"),
                ("would be drawn", f"{len(selection.kept):,}"),
                ("in lanes", f"{layout.lanes:,}"),
            ],
            indent="  ",
        )
        + "\n\n"
    )
    stdout.write(_preview(columns, body) + "\n")
    return 0


def run_types(arguments: argparse.Namespace, stdout: IO[str]) -> int:
    """List the feature kinds and their colours."""
    palette = Palette(arguments.palette, dict(arguments.colour))
    body = tabular.kind_rows(palette)
    if arguments.format == "tsv":
        stdout.write(tabular.render(tabular.KIND_COLUMNS, body))
        return 0
    if arguments.format == "json":
        columns = list(tabular.KIND_COLUMNS)
        stdout.write(
            json.dumps(
                {
                    "palette": palette.name,
                    "kinds": [dict(zip(columns, row)) for row in body],
                },
                indent=2,
            )
            + "\n"
        )
        return 0
    stdout.write(f"the {palette.name} palette, {len(body):,} kinds\n")
    stdout.write(
        prose.table(
            ["kind", "group", "colour"],
            [row[:3] for row in body],
            indent="  ",
        )
        + "\n"
    )
    return 0


def main(
    argv: Sequence[str] | None = None,
    stdout: IO[str] | None = None,
    stderr: IO[str] | None = None,
) -> int:
    """Run the program and return the code the shell should see.

    Errors are written to the handle given rather than raised, so that calling
    this from a test or another program behaves the way running it does.
    """
    out = stdout if stdout is not None else sys.stdout
    err = stderr if stderr is not None else sys.stderr
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        if arguments.command in ("circular", "linear"):
            return run_draw(arguments, out, arguments.command)
        if arguments.command == "features":
            return run_features(arguments, out)
        return run_types(arguments, out)
    except (
        UsageError,
        ReadError,
        DrawingError,
        FeatureError,
        ColourError,
        CompositionError,
    ) as error:
        err.write(f"plasmid-map: {error}\n")
        return 2
    except BrokenPipeError:
        return 0
