"""What colour each kind of feature is drawn.

Colour carries meaning on a plasmid map - coding sequences read as one thing,
origins and resistance markers as another - so the colours are grouped by what
a feature is for rather than spread evenly around a wheel. Kinds that answer
the same question share a hue and differ in lightness, which keeps a map
readable when it is printed in black and white and when it is read by someone
who cannot separate red from green.

There are two palettes. ``colour`` is the one above; ``mono`` is a set of greys
for a figure that will be printed without colour, and keeps the same grouping
so that the two versions of a map have the same structure.

Any kind not named here falls back through its group to a default, so a map of
a file using feature keys this module has never seen is drawn rather than
refused; the summary says which kinds took the default.
"""

from __future__ import annotations

import re

#: What a kind is for, which is what decides its hue.
GROUPS: dict[str, str] = {
    "CDS": "coding",
    "gene": "coding",
    "mat_peptide": "coding",
    "sig_peptide": "coding",
    "transit_peptide": "coding",
    "proprotein": "coding",
    "mRNA": "transcript",
    "precursor_RNA": "transcript",
    "exon": "transcript",
    "intron": "transcript",
    "5'UTR": "transcript",
    "3'UTR": "transcript",
    "ncRNA": "transcript",
    "rRNA": "transcript",
    "tRNA": "transcript",
    "tmRNA": "transcript",
    "misc_RNA": "transcript",
    "regulatory": "regulatory",
    "promoter": "regulatory",
    "enhancer": "regulatory",
    "RBS": "regulatory",
    "polyA_signal": "regulatory",
    "polyA_site": "regulatory",
    "TATA_signal": "regulatory",
    "minus_35_signal": "regulatory",
    "minus_10_signal": "regulatory",
    "GC_signal": "regulatory",
    "attenuator": "regulatory",
    "terminator": "terminator",
    "rep_origin": "origin",
    "oriT": "origin",
    "D-loop": "origin",
    "repeat_region": "repeat",
    "LTR": "repeat",
    "mobile_element": "repeat",
    "STS": "repeat",
    "primer_bind": "binding",
    "protein_bind": "binding",
    "misc_binding": "binding",
    "operator": "binding",
    "variation": "variant",
    "misc_difference": "variant",
    "old_sequence": "variant",
    "modified_base": "variant",
    "unsure": "variant",
    "stem_loop": "structure",
    "misc_structure": "structure",
    "source": "whole",
    "misc_feature": "other",
}

#: The colour of each group, and the fallback when a kind has no group.
COLOURS: dict[str, dict[str, str]] = {
    "colour": {
        "coding": "#3d6fb4",
        "transcript": "#5fa8c8",
        "regulatory": "#3f9a63",
        "terminator": "#c2553d",
        "origin": "#7f5aa8",
        "repeat": "#b07d3a",
        "binding": "#8c93a8",
        "variant": "#b8b3c4",
        "structure": "#55a39a",
        "whole": "#d8d8d8",
        "other": "#9aa0a6",
    },
    "mono": {
        "coding": "#3a3a3a",
        "transcript": "#6b6b6b",
        "regulatory": "#565656",
        "terminator": "#2a2a2a",
        "origin": "#808080",
        "repeat": "#9a9a9a",
        "binding": "#b4b4b4",
        "variant": "#cfcfcf",
        "structure": "#8a8a8a",
        "whole": "#e0e0e0",
        "other": "#c0c0c0",
    },
}

#: The palettes on offer.
PALETTES: tuple[str, ...] = tuple(COLOURS)

#: The group a kind with no group of its own is drawn as.
DEFAULT_GROUP: str = "other"


class ColourError(ValueError):
    """A colour that cannot be used."""


def group_of(kind: str) -> str:
    """Which group a feature kind belongs to.

    Unknown kinds fall into the default group rather than being refused, so a
    file using a feature key this module has not heard of still draws.
    """
    return GROUPS.get(kind, DEFAULT_GROUP)


def is_default(kind: str) -> bool:
    """Whether a kind is drawn in the default colour for want of its own."""
    return kind not in GROUPS


def check_colour(value: str) -> str:
    """A colour as written on the command line, if it can be used.

    Hexadecimal of three, four, six or eight digits, or one of the colour
    names every browser and every SVG renderer agrees on; anything else is
    refused here rather than producing a file that draws differently in two
    viewers.
    """
    text = value.strip()
    if re.fullmatch(r"#[0-9A-Fa-f]{3,8}", text) and len(text) in (4, 5, 7, 9):
        return text.lower()
    if re.fullmatch(r"[a-z]+", text) and text in NAMED:
        return text
    raise ColourError(
        f"not a colour: {value!r}; write it as #rrggbb or use a named colour"
    )


#: The colour names kept, which are the ones SVG 1.1 fixes exactly.
NAMED: frozenset[str] = frozenset(
    {
        "black",
        "white",
        "grey",
        "gray",
        "silver",
        "red",
        "maroon",
        "orange",
        "yellow",
        "olive",
        "green",
        "lime",
        "teal",
        "aqua",
        "cyan",
        "blue",
        "navy",
        "purple",
        "fuchsia",
        "magenta",
    }
)


class Palette:
    """The colour of every feature, with whatever overrides were asked for."""

    def __init__(
        self, name: str = "colour", overrides: dict[str, str] | None = None
    ) -> None:
        if name not in COLOURS:
            raise ColourError(
                f"not a palette: {name!r}; it is one of {', '.join(PALETTES)}"
            )
        self.name = name
        self.overrides = {
            kind: check_colour(value) for kind, value in (overrides or {}).items()
        }

    def colour(self, kind: str) -> str:
        """The colour one kind of feature is drawn in."""
        if kind in self.overrides:
            return self.overrides[kind]
        return COLOURS[self.name][group_of(kind)]

    def group_colour(self, group: str) -> str:
        """The colour of a whole group, for a legend."""
        return COLOURS[self.name].get(group, COLOURS[self.name][DEFAULT_GROUP])

    def rows(self) -> list[tuple[str, str, str]]:
        """Every kind this module knows, its group and its colour."""
        listed = sorted(GROUPS.items(), key=lambda pair: (pair[1], pair[0]))
        return [(kind, group, self.colour(kind)) for kind, group in listed]
