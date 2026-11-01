#!/usr/bin/env python3
"""Rebuild the example files, or check the ones committed here.

Three of the four examples are real records, fetched from the NCBI Nucleotide
database and written out unchanged. Each was chosen for something a plasmid
map has to get right:

* ``pbr322.gb`` (J01749.1) is the classic cloning vector: fifty features on
  4,361 bases, several of them one base long, several enclosing others.
* ``phix174.gb`` (NC_001422.1) has three genes written ``join(3981..5386,
  1..136)`` and the like, which cross the origin of a circular molecule. A map
  that draws those as two features, or gives them the length of the whole
  genome, is wrong, and these are what prove it is not.
* ``aav2.gb`` (NC_001401.2) is a linear genome with inverted terminal repeats
  at both ends and spliced coding sequences in the middle, which is what the
  linear track is for.

The fourth, ``cassette.tsv``, is not a record and is not fetched. It is a
feature table written by hand for an AAV transfer cassette, and it is here
because a table and a length are the other way into this program: no sequence,
no GenBank file, just what is where. It is committed rather than built, so it
is checked here but never rewritten.

Run with no arguments to fetch; run with ``--check`` to verify what is here
against the checksums below without touching the network.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent

#: Where the records come from. E-utilities serves the GenBank flat file
#: itself, which is what this program reads.
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

#: The accessions, pinned to a version so that a rebuild years from now
#: fetches the record these checksums were taken from rather than whatever
#: has replaced it.
RECORDS: dict[str, str] = {
    "pbr322.gb": "J01749.1",
    "phix174.gb": "NC_001422.1",
    "aav2.gb": "NC_001401.2",
}

#: The file that is written by hand rather than fetched.
COMPOSED: tuple[str, ...] = ("cassette.tsv",)

#: SHA-256 of each example as committed. A record edited by hand, or a
#: database that has changed one under us, fails --check rather than quietly
#: changing what the README's quoted output should say.
CHECKSUMS: dict[str, str] = {
    "pbr322.gb": "7660b0fe95a91d492b102192dc72924508a69fda762988aa0923b3d5cc2073f6",
    "phix174.gb": "0a8099e207633f8c3ba998bfe2a198518e1de356c7956691101b2512324b807c",
    "aav2.gb": "712b98e5f102572bcd3dbd7dcfb893eaa73350f6e47c6ab466ec30037c41cb8c",
    "cassette.tsv": "422e1963c8674a29542e2bfda45546bdccfa23bcc43df97a23234f3f711b9bac",
}


def digest(path: Path) -> str:
    """The SHA-256 of a file, as hexadecimal."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(accession: str) -> bytes:
    """One GenBank flat file from the NCBI Nucleotide database."""
    query = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": accession,
            "rettype": "gbwithparts",
            "retmode": "text",
            "tool": "plasmid-map-examples",
        }
    )
    with urllib.request.urlopen(f"{EUTILS}?{query}", timeout=120) as response:
        body = response.read()
    if not body.startswith(b"LOCUS"):
        raise SystemExit(
            f"{accession}: the service did not return a GenBank record"
        )
    return body


def check() -> int:
    """Verify every committed example against its recorded checksum."""
    failed = 0
    for name, expected in sorted(CHECKSUMS.items()):
        path = HERE / name
        if not path.exists():
            print(f"missing: {name}")
            failed += 1
            continue
        if not expected:
            print(f"no checksum recorded for {name}")
            failed += 1
            continue
        found = digest(path)
        if found != expected:
            print(f"changed: {name}\n  expected {expected}\n  found    {found}")
            failed += 1
    if failed:
        print(f"{failed} example(s) do not match what is recorded here")
        return 1
    print(f"{len(CHECKSUMS)} examples match their checksums")
    return 0


def rebuild() -> int:
    """Fetch the records again and report what the checksums should be."""
    for name, accession in sorted(RECORDS.items()):
        body = fetch(accession)
        (HERE / name).write_bytes(body)
        print(f"wrote {name} from {accession}")
    for name in sorted(CHECKSUMS):
        print(f'    "{name}": "{digest(HERE / name)}",')
    print(
        "Paste those lines into CHECKSUMS above if the records have "
        "legitimately changed, and re-run the test suite: the README quotes "
        "output measured from these files."
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    """Fetch the examples, or check the ones that are here."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify the committed examples against their checksums, offline",
    )
    arguments = parser.parse_args(argv)
    if arguments.check:
        return check()
    return rebuild()


if __name__ == "__main__":
    sys.exit(main())
