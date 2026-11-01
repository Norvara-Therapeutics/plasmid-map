"""Getting a construct out of a file, and what the file did not say."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from plasmid_map.reader import (
    ReadError,
    annotate,
    blank,
    detect,
    load,
    parse_fasta,
    read_text,
)

FASTA = ">cassette An example\nACGTACGTAC\nGTACGTACGT\n"


class TestDetection:
    def test_a_genbank_file_begins_with_a_locus_line(self, examples: Path) -> None:
        text = (examples / "pbr322.gb").read_text(encoding="utf-8")
        assert detect(text) == "genbank"

    def test_a_fasta_file_begins_with_a_greater_than(self) -> None:
        assert detect(FASTA) == "fasta"

    def test_anything_else_is_refused_with_what_was_expected(self) -> None:
        with pytest.raises(ReadError) as raised:
            detect("chr1\t10\t20\n")
        assert "LOCUS" in str(raised.value)


class TestFasta:
    def test_the_name_and_description_are_split_at_the_first_space(self) -> None:
        name, description, sequence = parse_fasta(FASTA)[0]
        assert name == "cassette"
        assert description == "An example"
        assert len(sequence) == 20

    def test_several_records_are_read(self) -> None:
        assert len(parse_fasta(FASTA + ">second\nACGT\n")) == 2

    def test_a_file_that_does_not_begin_with_a_heading_is_refused(self) -> None:
        with pytest.raises(ReadError):
            parse_fasta("ACGTACGT\n")

    def test_an_empty_file_is_refused(self) -> None:
        with pytest.raises(ReadError):
            parse_fasta("\n\n")


class TestLoading:
    def test_a_genbank_record_states_its_topology(self, examples: Path) -> None:
        loaded = load(str(examples / "pbr322.gb"))
        assert loaded.topology_stated
        assert loaded.first.is_circular

    def test_a_fasta_record_does_not(self, tmp_path: Path) -> None:
        path = tmp_path / "x.fasta"
        path.write_text(FASTA, encoding="utf-8")
        loaded = load(str(path))
        assert not loaded.topology_stated
        assert loaded.first.topology == "linear"

    def test_naming_the_topology_makes_it_stated(self, tmp_path: Path) -> None:
        path = tmp_path / "x.fasta"
        path.write_text(FASTA, encoding="utf-8")
        loaded = load(str(path), topology="circular")
        assert loaded.topology_stated
        assert loaded.first.is_circular

    def test_a_gzipped_file_is_read_whatever_it_is_called(
        self, tmp_path: Path, examples: Path
    ) -> None:
        body = (examples / "pbr322.gb").read_bytes()
        path = tmp_path / "pbr322.gb"
        path.write_bytes(gzip.compress(body))
        assert load(str(path)).first.length == 4361

    def test_a_file_that_is_not_there_is_refused_by_name(self) -> None:
        with pytest.raises(ReadError) as raised:
            load("/no/such/file.gb")
        assert "/no/such/file.gb" in str(raised.value)

    def test_a_directory_is_refused(self, tmp_path: Path) -> None:
        with pytest.raises(ReadError):
            load(str(tmp_path))

    def test_an_empty_file_is_refused(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.gb"
        path.write_text("", encoding="utf-8")
        with pytest.raises(ReadError):
            load(str(path))

    def test_a_format_that_does_not_exist_is_refused_by_name(
        self, examples: Path
    ) -> None:
        with pytest.raises(ReadError) as raised:
            load(str(examples / "pbr322.gb"), input_format="embl")
        assert "genbank" in str(raised.value)

    def test_further_records_are_counted_rather_than_lost(
        self, tmp_path: Path, examples: Path
    ) -> None:
        body = (examples / "pbr322.gb").read_text(encoding="utf-8")
        path = tmp_path / "two.gb"
        path.write_text(body + body, encoding="utf-8")
        loaded = load(str(path))
        assert len(loaded.constructs) == 2
        assert loaded.extra == 1

    def test_latin_one_is_read_when_utf_eight_will_not_decode(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "x.fasta"
        path.write_bytes(b">name\xe9\nACGT\n")
        assert load(str(path)).first.length == 4


class TestTopologyChange:
    def test_making_a_circular_record_linear_splits_what_crossed_the_origin(
        self, examples: Path
    ) -> None:
        """On a linear molecule it is two pieces, and saying otherwise would
        draw a feature running off the end."""
        loaded = load(str(examples / "phix174.gb"), topology="linear")
        construct = loaded.first
        assert not construct.is_circular
        assert construct.crossing() == ()
        gene = next(
            item
            for item in construct.features
            if item.kind == "gene" and item.start == 0 and len(item.parts) == 2
        )
        assert gene.length == 1542


class TestBlank:
    def test_a_length_alone_makes_a_construct(self) -> None:
        loaded = blank(2744, "linear", "cassette")
        assert loaded.first.length == 2744
        assert loaded.first.sequence is None
        assert not loaded.topology_stated

    def test_a_construct_needs_a_length(self) -> None:
        with pytest.raises(ReadError):
            blank(0)

    def test_a_topology_that_does_not_exist_is_refused(self) -> None:
        with pytest.raises(ReadError):
            blank(100, "round")


class TestAnnotate:
    def test_a_table_is_added_to_what_the_record_already_holds(
        self, examples: Path
    ) -> None:
        loaded = load(str(examples / "pbr322.gb"))
        before = len(loaded.first.features)
        after, added = annotate(
            loaded.first, typed=["insert:100..200:+:CDS"]
        )
        assert added == 1
        assert len(after.features) == before + 1

    def test_a_table_file_is_read(self, tmp_path: Path) -> None:
        path = tmp_path / "features.tsv"
        path.write_text("name\tstart\tend\nx\t10\t20\n", encoding="utf-8")
        construct = blank(100).first
        after, added = annotate(construct, table=str(path))
        assert added == 1
        assert after.features[0].name == "x"

    def test_nothing_to_add_leaves_the_construct_alone(self) -> None:
        construct = blank(100).first
        after, added = annotate(construct)
        assert added == 0
        assert after is construct

    def test_a_bad_table_is_refused_with_the_path_in_the_message(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "features.tsv"
        path.write_text("nothing\n", encoding="utf-8")
        with pytest.raises(ReadError) as raised:
            annotate(blank(100).first, table=str(path))
        assert str(path) in str(raised.value)

    def test_an_annotation_off_the_end_is_refused(self) -> None:
        with pytest.raises(ReadError):
            annotate(blank(100).first, typed=["x:90..200"])


class TestReadText:
    def test_a_missing_file_is_refused(self) -> None:
        with pytest.raises(ReadError):
            read_text("/no/such/file")

    def test_something_that_begins_like_gzip_but_is_not_is_refused(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "x.gz"
        path.write_bytes(b"\x1f\x8bnot really")
        with pytest.raises(ReadError) as raised:
            read_text(str(path))
        assert "gzipped" in str(raised.value)
