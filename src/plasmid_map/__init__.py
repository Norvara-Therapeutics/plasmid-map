"""Circular plasmid maps and linear feature tracks as standalone SVG.

Positions are held 0-based and half-open inside the program and written
1-based with both ends included everywhere a person or another program sees
them, which is what GenBank writes and what a sequence viewer shows. The
conversion happens in two places only - the readers on the way in, and
``rows`` on the way out - so nothing in between sees two conventions.

A feature that crosses the origin of a circular construct is one feature, held
as a single span whose end runs past the length of the construct. Its length
is the sum of its parts. Drawn as a circle it is one arc through twelve
o'clock; drawn as a line, which cuts the circle, it is the two pieces it
becomes, and the summary says which features that happened to.

Nothing is dropped quietly. A feature too short to see is drawn wider and
said to have been; a label with nowhere to go is moved along its side of the
map and joined to its feature by a leader line rather than left out; a feature
left out by a filter is counted and named. The picture cannot say any of that,
so the summary does, every time.

Every drawing is self-contained: no script, no web font, nothing fetched.
Nothing depends on the clock or on the order of a set, so the same construct
drawn twice gives the same bytes.
"""

from plasmid_map.composition import (
    CompositionError,
    Profile,
    Window,
    at_rich,
    gc_share,
    gc_skew,
    measure,
)
from plasmid_map.drawing import (
    Drawing,
    DrawingError,
    Selection,
    Style,
    groups_of,
    select,
)
from plasmid_map.feature import (
    STRANDS,
    TOPOLOGIES,
    Construct,
    Feature,
    FeatureError,
    Span,
    joined,
    location_text,
    minimum_width,
    overlapping,
    sorted_features,
)
from plasmid_map.genbank import GenBankError, Location, Locus, parse_location
from plasmid_map.genbank import parse as parse_genbank
from plasmid_map.labels import (
    Label,
    LabelError,
    Placed,
    Rows,
    place_circular,
    separate,
    shorten,
    stack_linear,
)
from plasmid_map.lanes import Layout, Placement, assign
from plasmid_map.palette import PALETTES, ColourError, Palette, group_of
from plasmid_map.reader import Loaded, ReadError, annotate, detect, load
from plasmid_map.ticks import RulerError, Tick, ruler, tick_step
from plasmid_map.tables import Annotation, TableError, parse_annotation, seqids

__version__ = "1.0.0rc2"

__all__ = [
    "Annotation",
    "ColourError",
    "CompositionError",
    "Construct",
    "Drawing",
    "DrawingError",
    "Feature",
    "FeatureError",
    "GenBankError",
    "Label",
    "LabelError",
    "Layout",
    "Loaded",
    "Location",
    "Locus",
    "PALETTES",
    "Palette",
    "Placed",
    "Placement",
    "Profile",
    "ReadError",
    "Rows",
    "RulerError",
    "STRANDS",
    "Selection",
    "Span",
    "Style",
    "TOPOLOGIES",
    "TableError",
    "Tick",
    "Window",
    "__version__",
    "annotate",
    "assign",
    "at_rich",
    "detect",
    "gc_share",
    "gc_skew",
    "group_of",
    "groups_of",
    "joined",
    "load",
    "location_text",
    "measure",
    "minimum_width",
    "overlapping",
    "parse_annotation",
    "parse_genbank",
    "parse_location",
    "place_circular",
    "ruler",
    "select",
    "separate",
    "seqids",
    "shorten",
    "sorted_features",
    "stack_linear",
    "tick_step",
]
