"""Feature tables and command-line annotations."""

from __future__ import annotations

import pytest

from plasmid_map.tables import (
    TableError,
    annotations,
    describe,
    looks_like_gff3,
    parse_annotation,
    read,
    read_gff3,
    read_table,
)

TSV = """name\tkind\tlocation\tnote
left ITR\trepeat_region\t1..145\tinverted terminal repeat
EGFP\tCDS\t1021..1740\t
primer\tprimer_bind\tcomplement(1760..1780)\treads backwards
"""

CSV = "name,type,start,end,strand\nbla,CDS,100,1000,-\ntet,CDS,1200,2000,+\n"

GFF = """##gff-version 3
p1\tplasmid-map\tCDS\t100\t1000\t.\t-\t0\tID=bla;Name=bla
p1\tplasmid-map\tpromoter\t50\t90\t.\t+\t.\tName=P1
"""


class TestDelimited:
    def test_a_table_counts_from_1_with_both_ends_in(self) -> None:
        features = read_table(TSV, 3000, False)
        itr = features[0]
        assert itr.parts[0].start == 0
        assert itr.parts[0].end == 145
        assert itr.length == 145

    def test_a_location_column_carries_a_strand(self) -> None:
        features = read_table(TSV, 3000, False)
        assert features[2].strand == -1

    def test_start_and_end_columns_are_read_too(self) -> None:
        features = read_table(CSV, 3000, False)
        assert [item.name for item in features] == ["bla", "tet"]
        assert features[0].strand == -1

    def test_the_delimiter_is_worked_out_from_the_header(self) -> None:
        assert len(read_table(CSV, 3000, False)) == 2
        assert len(read_table(TSV, 3000, False)) == 3

    def test_a_column_this_program_has_no_use_for_is_carried_through(self) -> None:
        features = read_table(TSV, 3000, False)
        assert features[0].qualifier("note") == "inverted terminal repeat"

    def test_a_row_with_no_kind_gets_the_default(self) -> None:
        features = read_table("name\tstart\tend\nx\t10\t20\n", 100, False)
        assert features[0].kind == "misc_feature"

    def test_a_table_with_no_position_columns_is_refused_by_name(self) -> None:
        with pytest.raises(TableError) as raised:
            read_table("name\tnote\nx\ty\n", 100, False)
        assert "location column" in str(raised.value)

    def test_a_table_with_no_delimiter_is_refused(self) -> None:
        with pytest.raises(TableError):
            read_table("nothing here\n", 100, False)

    def test_an_empty_table_is_refused(self) -> None:
        with pytest.raises(TableError):
            read_table("\n\n", 100, False)

    def test_a_position_that_is_not_a_number_is_refused_with_its_line(self) -> None:
        with pytest.raises(TableError) as raised:
            read_table("name\tstart\tend\nx\tfirst\t20\n", 100, False)
        assert "line 2" in str(raised.value)

    def test_counting_from_zero_is_refused(self) -> None:
        with pytest.raises(TableError) as raised:
            read_table("name\tstart\tend\nx\t0\t20\n", 100, False)
        assert "counts from 1" in str(raised.value)

    def test_an_end_before_a_start_is_refused(self) -> None:
        with pytest.raises(TableError):
            read_table("name\tstart\tend\nx\t50\t20\n", 100, False)

    def test_a_join_at_the_origin_becomes_one_feature(self) -> None:
        features = read_table(
            "name\tlocation\nwrap\tjoin(900..1000,1..100)\n", 1000, True
        )
        assert len(features[0].parts) == 1
        assert features[0].length == 201

    def test_the_two_ways_of_giving_a_position_say_different_things_about_strand(
        self,
    ) -> None:
        """A GenBank location without complement() is the plus strand; a start
        and an end column say nothing about strand at all."""
        located = read_table("name\tlocation\nx\t10..20\n", 100, False)
        columned = read_table("name\tstart\tend\nx\t10\t20\n", 100, False)
        assert located[0].strand == 1
        assert columned[0].strand == 0
        assert located[0].parts == columned[0].parts

    def test_a_strand_column_settles_it_either_way(self) -> None:
        both = read_table(
            "name\tlocation\tstrand\nx\t10..20\t-\ny\t30..40\t.\n", 100, False
        )
        assert [item.strand for item in both] == [-1, 0]


class TestGff3:
    def test_a_gff3_file_is_recognised(self) -> None:
        assert looks_like_gff3(GFF)
        assert not looks_like_gff3(CSV)

    def test_it_counts_from_1_with_both_ends_in(self) -> None:
        features = read_gff3(GFF, 3000, False)
        assert features[0].parts[0].start == 99
        assert features[0].parts[0].end == 1000

    def test_the_name_attribute_becomes_the_label(self) -> None:
        features = read_gff3(GFF, 3000, False)
        assert [item.name for item in features] == ["bla", "P1"]

    def test_the_third_column_becomes_the_kind(self) -> None:
        assert [item.kind for item in read_gff3(GFF, 3000, False)] == [
            "CDS",
            "promoter",
        ]

    def test_attributes_are_carried_through(self) -> None:
        features = read_gff3(GFF, 3000, False)
        assert features[0].qualifier("ID") == "bla"

    def test_comments_and_blank_lines_are_skipped(self) -> None:
        assert len(read_gff3("##gff-version 3\n\n# a comment\n", 100, False)) == 0

    def test_a_short_line_is_refused_with_its_line_number(self) -> None:
        with pytest.raises(TableError) as raised:
            read_gff3("##gff-version 3\np1\tx\tCDS\t1\t10\n", 100, False)
        assert "line 2" in str(raised.value)


class TestDispatch:
    def test_the_format_is_worked_out_when_it_is_not_given(self) -> None:
        assert len(read(GFF, 3000, False)) == 2
        assert len(read(TSV, 3000, False)) == 3

    def test_a_format_that_does_not_exist_is_refused_by_name(self) -> None:
        with pytest.raises(TableError) as raised:
            read(TSV, 100, False, kind="bed")
        assert "gff3" in str(raised.value)

    def test_every_format_is_described(self) -> None:
        described = dict(describe())
        assert set(described) == {"table", "gff3"}
        assert all(text for text in described.values())


class TestAnnotations:
    def test_a_simple_annotation_is_read(self) -> None:
        found = parse_annotation("EGFP:1021..1740")
        assert found.name == "EGFP"
        assert found.start == 1020
        assert found.end == 1740

    def test_a_strand_and_a_kind_may_follow(self) -> None:
        found = parse_annotation("EGFP:1021..1740:+:CDS")
        assert found.strand == 1
        assert found.kind == "CDS"

    def test_a_single_position_is_one_base(self) -> None:
        found = parse_annotation("cut:1200")
        assert found.end - found.start == 1

    def test_a_name_may_hold_a_colon(self) -> None:
        """The positions are found first, so everything before them is name."""
        found = parse_annotation("promoter:P1:100..200")
        assert found.name == "promoter:P1"

    def test_something_with_no_positions_is_refused_by_example(self) -> None:
        with pytest.raises(TableError) as raised:
            parse_annotation("just a name")
        assert "NAME:START..END" in str(raised.value)

    def test_an_annotation_with_no_name_is_refused(self) -> None:
        with pytest.raises(TableError):
            parse_annotation("100..200")

    def test_annotations_become_features(self) -> None:
        features = annotations(["a:1..10", "b:20..30:-"], 100, False)
        assert [item.name for item in features] == ["a", "b"]
        assert features[1].strand == -1

    def test_an_annotation_past_a_linear_end_is_refused(self) -> None:
        with pytest.raises(TableError) as raised:
            annotations(["a:90..200"], 100, False)
        assert "100" in str(raised.value)

    def test_an_annotation_past_a_circular_end_is_allowed(self) -> None:
        assert len(annotations(["a:90..120"], 100, True)) == 1


class TestExample:
    def test_the_committed_cassette_reads_as_nine_features(self, examples) -> None:
        text = (examples / "cassette.tsv").read_text(encoding="utf-8")
        features = read(text, 2744, False)
        assert len(features) == 9
        assert {item.kind for item in features} >= {"CDS", "promoter", "repeat_region"}
        assert any(item.strand == -1 for item in features)
