"""Reading GenBank flat files, including the parts real records get wrong."""

from __future__ import annotations

import pytest

from plasmid_map.feature import Construct
from plasmid_map.genbank import (
    GenBankError,
    parse,
    parse_location,
    parse_locus,
    parse_record,
    split_records,
)

RECORD = """LOCUS       TESTPLASMID             2000 bp    DNA     circular SYN 01-JAN-2026
DEFINITION  A test plasmid for reading, with a definition that
            runs onto a second line.
ACCESSION   X00001
VERSION     X00001.1
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            other sequences; artificial sequences.
FEATURES             Location/Qualifiers
     source          1..2000
                     /organism="synthetic construct"
     CDS             join(1900..2000,1..100)
                     /gene="wrap"
                     /product="a protein written across the origin"
     CDS             complement(500..800)
                     /gene="rev"
     misc_feature    1200
                     /note="one base"
     misc_feature    1300^1301
                     /note="between two bases"
     regulatory      <1500..>1600
                     /regulatory_class="promoter"
                     /note="an end that runs off the record"
ORIGIN
{body}//
"""

#: A sequence as long as the LOCUS line of the sample record says it is.
BASES: str = ("acgt" * 500)[:2000]


def origin(sequence: str) -> str:
    """An ORIGIN block holding a sequence, written as GenBank writes one."""
    lines = []
    for start in range(0, len(sequence), 60):
        blocks = [
            sequence[position : position + 10]
            for position in range(start, min(start + 60, len(sequence)), 10)
        ]
        lines.append(f"{start + 1:>9} " + " ".join(blocks))
    return "\n".join(lines) + ("\n" if lines else "")


def sample(sequence: str = BASES) -> str:
    """The sample record carrying a sequence of a given length."""
    return RECORD.format(body=origin(sequence))


class TestLocus:
    def test_the_fields_are_read_from_the_left(self) -> None:
        locus = parse_locus(
            "LOCUS       NC_001422               5386 bp ss-DNA     "
            "circular PHG 11-JAN-2023"
        )
        assert locus.name == "NC_001422"
        assert locus.length == 5386
        assert locus.molecule == "ss-DNA"
        assert locus.is_circular

    def test_a_line_that_stops_early_has_no_topology_rather_than_a_wrong_one(
        self,
    ) -> None:
        """Reading right to left would make 'DNA' the topology."""
        locus = parse_locus("LOCUS       SHORT                    100 bp    DNA")
        assert locus.topology == "linear"
        assert locus.molecule == "DNA"

    def test_a_line_with_no_length_is_refused(self) -> None:
        with pytest.raises(GenBankError):
            parse_locus("LOCUS       SHORT    bp    DNA  circular")

    def test_a_protein_record_is_refused_by_name(self) -> None:
        with pytest.raises(GenBankError) as raised:
            parse_locus("LOCUS       NP_040703                427 aa            linear")
        assert "protein" in str(raised.value)

    def test_something_that_is_not_a_locus_line_is_refused(self) -> None:
        with pytest.raises(GenBankError):
            parse_locus("DEFINITION  not a locus line")


class TestLocations:
    def test_a_range_is_converted_once(self) -> None:
        where = parse_location("86..1276")
        assert where.parts[0].start == 85
        assert where.parts[0].end == 1276

    def test_a_single_base_is_one_base_long(self) -> None:
        assert parse_location("426").parts[0].length == 1

    def test_a_site_between_two_bases_has_no_length(self) -> None:
        assert parse_location("102^103").parts[0].length == 0

    def test_complement_flips_the_strand(self) -> None:
        assert parse_location("complement(3293..4153)").strand == -1

    def test_a_join_keeps_its_parts(self) -> None:
        where = parse_location("join(321..1906,2228..2252)")
        assert [(part.start, part.end) for part in where.parts] == [
            (320, 1906),
            (2227, 2252),
        ]

    def test_order_is_read_like_join(self) -> None:
        assert len(parse_location("order(1..10,20..30)").parts) == 2

    def test_a_complement_of_a_join_is_one_feature_on_the_minus_strand(self) -> None:
        where = parse_location("complement(join(1..10,20..30))")
        assert where.strand == -1
        assert len(where.parts) == 2

    def test_partial_ends_are_recorded(self) -> None:
        assert parse_location("<1..>5386").partial == (True, True)

    def test_a_complement_swaps_which_end_is_partial(self) -> None:
        assert parse_location("complement(<1..300)").partial == (False, True)

    def test_a_location_in_another_record_is_refused_by_name(self) -> None:
        with pytest.raises(GenBankError) as raised:
            parse_location("J00194.1:100..202")
        assert "another record" in str(raised.value)

    def test_several_possible_positions_are_refused(self) -> None:
        with pytest.raises(GenBankError) as raised:
            parse_location("one-of(1888,1901)")
        assert "several" in str(raised.value)

    def test_a_location_that_ends_before_it_starts_is_refused(self) -> None:
        with pytest.raises(GenBankError):
            parse_location("500..400")

    def test_an_empty_location_is_refused(self) -> None:
        with pytest.raises(GenBankError):
            parse_location("   ")


@pytest.fixture(scope="module")
def record() -> Construct:
    """The sample record above, read once."""
    return parse(sample())[0]


class TestRecord:
    def test_the_definition_is_folded_onto_one_line(self, record: Construct) -> None:
        assert record.description.endswith("runs onto a second line")

    def test_the_accession_is_read(self, record: Construct) -> None:
        assert record.accession == "X00001"

    def test_the_source_is_read(self, record: Construct) -> None:
        assert record.source.startswith("synthetic construct")

    def test_every_feature_is_read(self, record: Construct) -> None:
        assert len(record.features) == 6

    def test_a_join_at_the_origin_is_one_feature(self, record: Construct) -> None:
        wrapped = next(item for item in record.features if item.name == "wrap")
        assert len(wrapped.parts) == 1
        assert wrapped.length == 201
        assert wrapped.crosses_origin(2000)

    def test_a_feature_takes_its_name_from_the_gene_qualifier(
        self, record: Construct
    ) -> None:
        assert any(item.name == "rev" for item in record.features)

    def test_a_feature_with_only_a_note_takes_the_note(
        self, record: Construct
    ) -> None:
        assert any(item.name == "one base" for item in record.features)

    def test_a_source_feature_takes_the_organism_not_the_note(self) -> None:
        """A /note on a source feature is usually a paragraph, not a label."""
        record = parse(sample())[0]
        source = next(item for item in record.features if item.kind == "source")
        assert source.name == "synthetic construct"

    def test_a_wrapped_qualifier_is_rejoined(self, record: Construct) -> None:
        wrapped = next(item for item in record.features if item.name == "wrap")
        assert wrapped.qualifier("product") == "a protein written across the origin"

    def test_partial_ends_survive_into_the_feature(self, record: Construct) -> None:
        partial = next(item for item in record.features if item.kind == "regulatory")
        assert partial.is_partial


class TestSequence:
    def test_a_declared_length_the_origin_block_does_not_hold_is_refused(
        self,
    ) -> None:
        """Naming both figures, rather than quietly believing one of them."""
        with pytest.raises(GenBankError) as raised:
            parse(sample("acgt"))
        assert "2,000" in str(raised.value) and "4" in str(raised.value)

    def test_a_record_with_no_origin_block_has_no_sequence(self) -> None:
        record = parse(sample("").replace("ORIGIN\n", ""))[0]
        assert record.sequence is None

    def test_the_sequence_is_read_without_its_position_numbers(self) -> None:
        record = parse(sample())[0]
        assert record.sequence is not None
        assert set(record.sequence) <= set("ACGT")


class TestFiles:
    def test_several_records_are_split_on_the_terminator(self) -> None:
        text = sample() * 3
        assert len(split_records(text)) == 3
        assert len(parse(text)) == 3

    def test_a_file_with_no_locus_line_is_refused(self) -> None:
        with pytest.raises(GenBankError):
            parse_record(["DEFINITION  nothing here"])

    def test_an_empty_file_is_refused(self) -> None:
        with pytest.raises(GenBankError):
            parse("")


class TestRealRecords:
    def test_pbr322_reads_as_the_record_describes_it(self, pbr322) -> None:
        assert pbr322.length == 4361
        assert pbr322.is_circular
        assert len(pbr322.features) == 50
        assert pbr322.sequence is not None

    def test_phix174_has_three_genes_across_the_origin(self, phix174) -> None:
        crossing = phix174.crossing()
        assert len({item.name for item in crossing}) == 3
        assert len(crossing) == 6

    def test_the_phix174_a_gene_is_not_the_whole_genome(self, phix174) -> None:
        """join(3981..5386,1..136) is 1,542 bases, not 5,386."""
        gene = next(
            item
            for item in phix174.crossing()
            if item.kind == "gene" and item.start == 3980
        )
        assert gene.length == 1542

    def test_aav2_is_linear_and_has_spliced_coding_sequences(self, aav2) -> None:
        assert not aav2.is_circular
        spliced = [item for item in aav2.features if len(item.parts) > 1]
        assert len(spliced) == 2
        assert all(item.kind == "CDS" for item in spliced)
