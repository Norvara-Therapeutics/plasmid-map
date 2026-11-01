"""Fixtures shared by the test suite."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from plasmid_map.feature import Construct, Feature, Span  # noqa: E402


@pytest.fixture(scope="session")
def examples() -> Path:
    """The directory holding the example files."""
    return ROOT / "examples"


@pytest.fixture(scope="session")
def pbr322(examples: Path) -> Construct:
    """pBR322: fifty real features, several of them one base long."""
    from plasmid_map.genbank import parse

    return parse((examples / "pbr322.gb").read_text(encoding="utf-8"))[0]


@pytest.fixture(scope="session")
def phix174(examples: Path) -> Construct:
    """phiX174: three genes written across the origin."""
    from plasmid_map.genbank import parse

    return parse((examples / "phix174.gb").read_text(encoding="utf-8"))[0]


@pytest.fixture(scope="session")
def aav2(examples: Path) -> Construct:
    """AAV2: a linear genome with repeats at both ends and spliced genes."""
    from plasmid_map.genbank import parse

    return parse((examples / "aav2.gb").read_text(encoding="utf-8"))[0]


def feature(
    kind: str = "CDS",
    name: str = "",
    start: int = 0,
    end: int = 100,
    strand: int = 0,
    parts: tuple[tuple[int, int], ...] | None = None,
) -> Feature:
    """One feature, written the short way tests want to write it."""
    spans = (
        tuple(Span(first, last) for first, last in parts)
        if parts
        else (Span(start, end),)
    )
    return Feature(kind=kind, name=name or kind, parts=spans, strand=strand)


def construct(
    length: int = 1000,
    topology: str = "circular",
    features: tuple[Feature, ...] = (),
    sequence: str | None = None,
    name: str = "test",
) -> Construct:
    """A construct, written the short way tests want to write it."""
    return Construct(
        name=name,
        length=length,
        topology=topology,
        features=features,
        sequence=sequence,
    )


def spread(count: int = 24, length: int = 1000, step: int = 7) -> tuple[Feature, ...]:
    """Features that overlap each other in the awkward ways real ones do.

    Deliberately over-represented: features starting at the same base, a
    feature enclosing several others, one-base features a few bases apart, and
    one crossing the origin.
    """
    made: list[Feature] = []
    for number in range(count):
        start = (number * step) % length
        made.append(
            feature(
                kind="CDS" if number % 3 else "misc_feature",
                name=f"f{number}",
                start=start,
                end=start + (1 if number % 4 == 0 else 40),
                strand=(1, -1, 0)[number % 3],
            )
        )
    made.append(feature(kind="gene", name="whole", start=0, end=length - 1))
    made.append(
        feature(kind="CDS", name="wrapped", start=length - 30, end=length + 25)
    )
    return tuple(made)


def sequence(length: int = 1000, seed: int = 7) -> str:
    """A reproducible sequence with its GC varying along it."""
    import random

    generator = random.Random(seed)
    letters: list[str] = []
    for position in range(length):
        weights = (3, 1, 1, 3) if position % 400 < 200 else (1, 3, 3, 1)
        letters.append(generator.choices("ACGT", weights=weights)[0])
    return "".join(letters)
