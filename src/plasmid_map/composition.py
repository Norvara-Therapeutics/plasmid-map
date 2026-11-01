"""GC content and GC skew, measured in windows along the sequence.

On a circular construct the windows run round the origin rather than stopping
at it, so the ring has no seam and no window is short. On a linear one they
stop, and the last window is whatever is left. Which of the two happened is
reported, because a GC track with a seam at the origin is the usual sign that
a circular sequence was measured as a linear one.

GC content is the share of A, C, G and T that is G or C; anything else - N, an
IUPAC ambiguity code, a gap - is left out of both sides of the division rather
than counted as not-GC, which would pull a window of mostly N towards zero and
make it look AT-rich. How much of the sequence was left out is reported.

GC skew is ``(G - C) / (G + C)``, which has no value in a window holding
neither G nor C; that window is reported as having none rather than as zero,
since no bases to divide and an even split are different facts.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

#: The bases a GC share is taken over. Anything else is left out of both the
#: numerator and the denominator.
COUNTED: frozenset[str] = frozenset("ACGTU")

#: The bases that count towards the numerator.
STRONG: frozenset[str] = frozenset("CG")

#: The most windows a profile is built from, so that an SVG stays a size a
#: browser and a figure editor will both open.
WINDOW_LIMIT: int = 2000


class CompositionError(ValueError):
    """A profile that cannot be measured as asked."""


@dataclass(frozen=True)
class Window:
    """One window of the sequence and what it holds.

    ``start`` and ``end`` are 0-based and half-open; on a circular construct
    ``end`` may run past the length, which is how a window crossing the origin
    is held, the same way an origin-crossing feature is.
    """

    start: int
    end: int
    gc: float | None
    skew: float | None
    counted: int
    skipped: int

    @property
    def middle(self) -> float:
        """The position the window is drawn at."""
        return (self.start + self.end) / 2


@dataclass(frozen=True)
class Profile:
    """A whole sequence measured in windows."""

    windows: tuple[Window, ...]
    window: int
    step: int
    mean: float | None
    skipped: int
    length: int
    circular: bool

    @property
    def lowest(self) -> float | None:
        """The lowest GC share any window reached."""
        shares = [item.gc for item in self.windows if item.gc is not None]
        return min(shares) if shares else None

    @property
    def highest(self) -> float | None:
        """The highest GC share any window reached."""
        shares = [item.gc for item in self.windows if item.gc is not None]
        return max(shares) if shares else None

    @property
    def skew_limit(self) -> float:
        """The largest skew reached either way, for scaling a skew track."""
        values = [abs(item.skew) for item in self.windows if item.skew is not None]
        return max(values) if values else 0.0


def gc_share(sequence: str) -> tuple[float | None, int, int]:
    """The GC share of a stretch of sequence, and what was counted.

    Returns the share, how many bases went into it, and how many were skipped
    for not being one of A, C, G, T or U.
    """
    counted = 0
    strong = 0
    skipped = 0
    for base in sequence.upper():
        if base in COUNTED:
            counted += 1
            if base in STRONG:
                strong += 1
        else:
            skipped += 1
    if not counted:
        return None, 0, skipped
    return strong / counted, counted, skipped


def gc_skew(sequence: str) -> float | None:
    """``(G - C) / (G + C)``, or None where there is neither G nor C."""
    upper = sequence.upper()
    guanine = upper.count("G")
    cytosine = upper.count("C")
    total = guanine + cytosine
    if not total:
        return None
    return (guanine - cytosine) / total


def default_window(length: int) -> int:
    """A window size that gives a readable track at this length.

    A fortieth of the construct, never finer than fifty bases: fine enough to
    show the AT-rich stretch an origin of replication sits in, coarse enough
    that the ring is a shape rather than a hedge.
    """
    return max(50, min(length, length // 40))


def default_step(window: int, length: int) -> int:
    """How far apart the windows are placed, given their size."""
    step = max(1, window // 5)
    while length // step > WINDOW_LIMIT:
        step *= 2
    return step


def _slice(sequence: str, start: int, end: int, length: int) -> str:
    """The bases of a window, taken round the origin when it crosses one."""
    if end <= length:
        return sequence[start:end]
    return sequence[start:length] + sequence[: end - length]


def measure(
    sequence: str,
    circular: bool,
    window: int | None = None,
    step: int | None = None,
) -> Profile:
    """GC content and skew along a sequence, window by window.

    On a circular sequence every window is the full width and the last one
    meets the first, so the track closes. On a linear one the windows stop at
    the end and the final window is short, which is reported rather than
    padded.
    """
    length = len(sequence)
    if not length:
        raise CompositionError("there is no sequence to measure")
    size = window if window is not None else default_window(length)
    if size <= 0:
        raise CompositionError(f"a window has to hold bases, not {size}")
    size = min(size, length)
    stride = step if step is not None else default_step(size, length)
    if stride <= 0:
        raise CompositionError(f"a step has to move along, not {stride}")
    built: list[Window] = []
    skipped_total = 0
    start = 0
    last = False
    while start < length:
        end = start + size
        if not circular and end >= length:
            end = length
            last = True
        piece = _slice(sequence, start, end, length)
        share, _counted, skipped = gc_share(piece)
        built.append(
            Window(
                start=start,
                end=end,
                gc=share,
                skew=gc_skew(piece),
                counted=_counted,
                skipped=skipped,
            )
        )
        skipped_total += skipped
        if last:
            break
        start += stride
    whole, _counted, whole_skipped = gc_share(sequence)
    return Profile(
        windows=tuple(built),
        window=size,
        step=stride,
        mean=whole,
        skipped=whole_skipped,
        length=length,
        circular=circular,
    )


def at_rich(profile: Profile, count: int = 1) -> list[Window]:
    """The windows with the least GC, lowest first.

    Useful enough to be worth naming: the AT-rich stretch of a plasmid is
    usually where its origin of replication is, so a map drawn without
    annotation still points at one.
    """
    measured = [item for item in profile.windows if item.gc is not None]
    return sorted(measured, key=lambda item: (item.gc, item.start))[:count]


def shares(windows: Sequence[Window]) -> list[float]:
    """The GC shares that exist, for a summary that has to average them."""
    return [item.gc for item in windows if item.gc is not None]
