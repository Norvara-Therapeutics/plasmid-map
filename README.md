# plasmid-map

Draw circular plasmid maps and linear feature tracks as standalone SVG, from
GenBank files, feature tables, or a few annotations typed on the command line.
No third-party dependencies, no compiled extensions, no account, no upload.

Every drawing is one file with nothing outside it: no script, no web font,
nothing fetched. Every feature is a shape with a class, so a map opens in
Illustrator or Inkscape as something to edit rather than a picture of
something. The same input draws the same bytes, so two maps can be compared
with `diff`.

## Install

```bash
pip install --pre norvara-plasmid-map
```

Python 3.10 or newer. Nothing else.

## Quickstart

A map of pBR322, restricted to the features worth putting on a figure, with a
GC track inside the ruler. The drawing goes to `--svg`; what goes to the
terminal is an account of what was drawn.

```bash
plasmid-map circular examples/pbr322.gb --kind CDS,rep_origin,regulatory,repeat_region --gc --width 900 --svg pbr322.svg
```

```
examples/pbr322.gb (genbank, circular, 4,361 bp) drawn as a circular map
  features in the record     50
  drawn                      14 in 3 lanes
  not drawn                  36
  crossing the origin        0
  drawn wider than they are  10, the narrowest drawn as 9 bp
  labels                     14 labels, 8 moved apart, furthest 9 px
  ruler                      a tick every 200 bp
  canvas                     900 x 536 px
  GC content                 53.7% overall, 26.6% to 68.8% across the windows
  windows                    208 windows of 109 bp every 21 bp, round the origin

  kind           drawn  colour
  -------------  -----  -------
  regulatory         6  #3f9a63
  repeat_region      4  #b07d3a
  CDS                3  #3d6fb4
  rep_origin         1  #7f5aa8

  note: Left out by the kind filter: gene (2), mat_peptide (1), misc_binding (15), misc_difference (3), misc_feature (5), mobile_element (1), old_sequence (6), sig_peptide (1), source (2).

  note: 10 features are shorter than 9 bp and are drawn that wide so they can be seen, so their width on the map is not their length.

  note: 8 labels were moved along their side of the map to keep them apart, the furthest by 9 px. Each is joined to what it names by a leader line; none was dropped.
```

Three of phiX174's genes cross the origin. The largest is written
`join(3981..5386,1..136)` and is 1,542 bases long - the sum of its parts, not
the 5,386 that subtracting its first coordinate from its last would give, and
not two genes. Each of the three is drawn as one arc through twelve o'clock.

```bash
plasmid-map circular examples/phix174.gb --kind gene
```

```
examples/phix174.gb (genbank, circular, 5,386 bp) drawn as a circular map
  features in the record  32
  drawn                   11 in 4 lanes
  not drawn               21
  crossing the origin     3
  labels                  11 labels, 4 moved apart, furthest 4 px
  ruler                   a tick every 500 bp
  canvas                  760 x 616 px

  kind  drawn  colour
  ----  -----  -------
  gene     11  #3d6fb4

  note: Left out by the kind filter: CDS (11), mRNA (2), misc_feature (2), rep_origin (1), source (1), variation (4).

  note: 3 features cross the origin and are each drawn as one arc through it, not two: phiX174p01, phiX174p02, phiX174p03. Their lengths are the sum of their parts.

  note: 4 labels were moved along their side of the map to keep them apart, the furthest by 4 px. Each is joined to what it names by a leader line; none was dropped.
```

Drawing the same record as a line cuts the circle at the origin, so those
three genes become the two pieces they are, one at each end of the track. The
summary says so rather than leaving it to be noticed.

```bash
plasmid-map linear examples/phix174.gb --kind gene
```

```
examples/phix174.gb (genbank, circular, 5,386 bp) drawn as a linear map
  features in the record  32
  drawn                   11 in 4 rows
  not drawn               21
  crossing the origin     3
  labels                  11 labels in 4 rows
  ruler                   a tick every 500 bp
  canvas                  760 x 228 px

  kind  drawn  colour
  ----  -----  -------
  gene     11  #3d6fb4

  note: Left out by the kind filter: CDS (11), mRNA (2), misc_feature (2), rep_origin (1), source (1), variation (4).

  note: 3 features cross the origin, and a linear track cuts there, so each is drawn as the two pieces it becomes, at either end: phiX174p01, phiX174p02, phiX174p03.

  note: The labels needed 4 rows to keep them apart; each is joined to what it names by a leader line.
```

There need not be a sequence at all. A cassette being designed is a name, a
length and a list of what is where, and `--length` with a feature table draws
it.

```bash
plasmid-map linear --length 2744 --features examples/cassette.tsv --title "AAV-CMV-EGFP cassette"
```

```
the command line (length, linear, 2,744 bp) drawn as a linear map
  features in the record  9
  drawn                   9 in 2 rows
  not drawn               0
  labels                  9 labels in 2 rows
  ruler                   a tick every 200 bp
  canvas                  760 x 161 px

  kind           drawn  colour
  -------------  -----  -------
  repeat_region      2  #b07d3a
  CDS                1  #3d6fb4
  enhancer           1  #3f9a63
  intron             1  #5fa8c8
  misc_feature       1  #9aa0a6
  polyA_signal       1  #3f9a63
  primer_bind        1  #8c93a8
  promoter           1  #3f9a63

  note: The input does not record whether the ends are joined, so it was drawn as linear. Say which it is with --topology.

  note: The labels needed 2 rows to keep them apart; each is joined to what it names by a leader line.

  note: 9 features were added to the record before drawing, from a feature table or the command line.
```

## The commands

| Command | What it does |
| --- | --- |
| `circular` | draws the construct as a circle, features on concentric rings |
| `linear` | draws it as a horizontal track, features in stacked rows |
| `features` | lists what would be drawn, after the filters, without drawing it |
| `types` | lists the feature kinds known and the colour each is given |

`circular` and `linear` take the same options, except `--start`, which turns a
circle and has nothing to turn on a line, and so is not offered there.

## Options

| Option | Effect |
| --- | --- |
| `--length BASES` | draw a feature table or `--annotate` with no sequence file |
| `--features PATH`, `--features-format`, `--features-seqid` | add features from a delimited table or GFF3 |
| `--annotate NAME:START..END[:STRAND[:KIND]]` | add one feature; repeatable |
| `--topology`, `--input-format` | say what the file does not say |
| `--record NAME-OR-NUMBER` | which record of a file holding several to draw |
| `--kind TYPE[,TYPE]`, `--hide`, `--minimum-length` | choose what is drawn |
| `--qualifiers KEY[,KEY]` | carry qualifiers across as extra table columns |
| `--start BASE` | which base sits at twelve o'clock (circular only) |
| `--gc`, `--skew`, `--gc-window`, `--gc-step` | the composition tracks |
| `--width`, `--track`, `--labels`, `--label-width` | the size and the labelling |
| `--palette`, `--colour KIND=COLOUR`, `--legend`, `--title`, `--no-ruler` | the look |
| `--svg PATH`, `--report PATH`, `-o PATH`, `-f` | where the results go |

`--svg` writes the drawing, `--report` a standalone HTML page holding the map,
the summary and the tables behind it, and `-o`/`--output` the drawn features
as a tab-separated table. `-f`/`--format` chooses what goes to standard
output: `table`, `tsv`, `json`, `svg` or `none`. Any input may be `-` for
standard input, and a gzipped file is read as one whatever it is called. The
`features` command takes `--shape` to say which map's lanes to work out.

Two outputs cannot share one stream, so `--svg -` alongside anything other
than `-f none` is refused rather than producing an SVG and a JSON document
run together. An option that would have no effect as things stand - a
`--gc-window` with no track to size, a `--features-format` with no table to
read - is refused for the same reason.

## Reading the output

**Nothing is dropped quietly.** A feature too short to see is drawn wider and
the summary says how many and how wide. A label with nowhere to go is moved
along its side of the map and joined to its feature by a leader line, never
left out, and the summary says how many moved and how far. A feature left out
by a filter is counted and named. The picture cannot say any of that, which is
why the summary always does.

**Labels do not overlap.** On a circle they are sorted by where they would
like to be, runs that collide are merged into a block, and each block is
centred on the average of its members' wishes, which is the arrangement that
moves them least and keeps them in the order their features are in. On a line
they are stacked into rows above the track. How wide a label will be is
*estimated*, from a fixed average glyph width, since there is no font here to
measure with; the same estimate places the labels and checks that they do not
collide, so the two always agree with each other.

**Coordinates count from 1 and include both ends**, in every table, every
tooltip and every message, which is what GenBank writes and what a sequence
viewer shows. A feature crossing the origin has an `end` smaller than its
`start`, as it does in a GenBank file, a `crosses_origin` column that says so,
and a `length` that is the sum of its parts.

**The map fits the width you ask for.** The labels are measured first, so the
circle is made as large as what is left allows rather than the canvas growing
by however long the longest gene name happened to be. When a construct needs
more lanes than will fit, the bands are thinned rather than reaching the
centre, and the summary says they were.

**Each table in the HTML report downloads** as the tab-separated file `-o`
would have written, carried as a data URI; both come from the same serialiser,
and a test asserts they are identical byte for byte.

## What it does not do

It does not edit sequence, simulate a cloning step, find restriction sites or
open reading frames, or align anything; it draws what it is told is there.
It reads GenBank and FASTA, and feature tables as delimited text or GFF3, but
not SnapGene, Vector NTI, ApE or EMBL files. It refuses a GenBank location
that points at another record or names several possible positions, rather than
guessing which one to draw. It draws SVG only: there is no PNG or PDF writer,
because every renderer and every figure editor already has one.

## Development

```bash
pip install -e ".[dev]"
python -m pytest
python -m pyflakes src/plasmid_map tests scripts examples
python scripts/check_wheel.py
```

`examples/fetch_examples.py` rebuilds the three GenBank examples from the NCBI
Nucleotide database; `--check` verifies the committed ones against their
recorded checksums without touching the network. `examples/cassette.tsv` is
not a record: it is a feature table written by hand, and is checked but never
rewritten.

## License

PolyForm Noncommercial License 1.0.0. Free for noncommercial use; commercial
use requires a separate licence from Norvara Therapeutics.
