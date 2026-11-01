"""The three real records, and the answers they are known to have.

These are the tests that would catch a change of convention anywhere in the
program: the coordinates, the lengths and the counts here come from the
records themselves, and every one of them is a number a reader could check
against the GenBank entry.
"""

from __future__ import annotations

import pytest

from plasmid_map.circular import draw
from plasmid_map.composition import measure
from plasmid_map.drawing import Style, select
from plasmid_map.linear import draw as draw_linear
from plasmid_map.rows import FEATURE_COLUMNS, feature_rows


class TestPbr322:
    def test_it_is_the_plasmid_it_says_it_is(self, pbr322) -> None:
        assert pbr322.length == 4361
        assert pbr322.is_circular
        assert pbr322.accession == "J01749"

    def test_the_tetracycline_resistance_gene_is_where_the_record_puts_it(
        self, pbr322
    ) -> None:
        tet = next(
            item
            for item in pbr322.features
            if item.kind == "CDS" and item.name == "tet"
        )
        assert (tet.start + 1, tet.end) == (86, 1276)
        assert tet.length == 1191
        assert tet.strand == 1

    def test_the_beta_lactamase_gene_reads_the_other_way(self, pbr322) -> None:
        bla = next(
            item
            for item in pbr322.features
            if item.kind == "CDS" and item.name == "bla"
        )
        assert (bla.start + 1, bla.end) == (3293, 4153)
        assert bla.strand == -1

    def test_the_origin_of_replication_is_one_base_in_the_record(
        self, pbr322
    ) -> None:
        origin = next(item for item in pbr322.features if item.kind == "rep_origin")
        assert origin.length == 1
        assert origin.start + 1 == 2535

    def test_its_gc_content_is_what_it_is_known_for(self, pbr322) -> None:
        profile = measure(pbr322.sequence, True)
        assert profile.mean == pytest.approx(0.538, abs=0.002)

    def test_a_one_base_feature_is_drawn_wide_enough_to_see(self, pbr322) -> None:
        drawing = draw(pbr322, select(pbr322), Style(width=900))
        widened = {item.feature.kind for item in drawing.layout.widened}
        assert "rep_origin" in widened

    def test_the_whole_map_draws_without_a_label_collision(self, pbr322) -> None:
        drawing = draw(pbr322, select(pbr322), Style(width=900))
        assert drawing.drawn == 48
        assert drawing.label_collisions == []


class TestPhix174:
    def test_it_is_the_genome_it_says_it_is(self, phix174) -> None:
        assert phix174.length == 5386
        assert phix174.is_circular

    def test_gene_a_is_one_feature_across_the_origin(self, phix174) -> None:
        """Written join(3981..5386,1..136) in the record."""
        gene = next(
            item
            for item in phix174.features
            if item.kind == "gene" and item.start == 3980
        )
        assert len(gene.parts) == 1
        assert gene.length == 1542
        assert gene.crosses_origin(5386)

    def test_its_table_row_reads_the_way_the_record_does(self, phix174) -> None:
        drawing = draw(
            phix174, select(phix174, kinds=["gene"]), Style(width=900)
        )
        rows = feature_rows(drawing.layout.placements, phix174.length)
        row = next(row for row in rows if row[2] == 3981)
        assert row[FEATURE_COLUMNS.index("end")] == 136
        assert row[FEATURE_COLUMNS.index("length")] == 1542
        assert row[FEATURE_COLUMNS.index("crosses_origin")] is True
        assert row[FEATURE_COLUMNS.index("location")] == "join(3981..5386,1..136)"

    def test_three_genes_cross_the_origin_and_three_coding_sequences_do(
        self, phix174
    ) -> None:
        crossing = phix174.crossing()
        assert sum(1 for item in crossing if item.kind == "gene") == 3
        assert sum(1 for item in crossing if item.kind == "CDS") == 3

    def test_drawn_as_a_circle_each_is_a_single_arc(self, phix174) -> None:
        drawing = draw(
            phix174, select(phix174, kinds=["gene"]), Style(width=900)
        )
        for item in drawing.layout.placements:
            assert len(item.spans) == 1

    def test_drawn_as_a_line_each_becomes_two_blocks(self, phix174) -> None:
        drawing = draw_linear(
            phix174, select(phix174, kinds=["gene"]), Style(width=900)
        )
        assert len(drawing.split) == 3
        assert drawing.svg.count('class="feature"') == len(
            drawing.layout.placements
        ) + 3


class TestAav2:
    def test_it_is_the_genome_it_says_it_is(self, aav2) -> None:
        assert aav2.length == 4679
        assert not aav2.is_circular

    def test_the_inverted_terminal_repeats_sit_at_both_ends(self, aav2) -> None:
        repeats = [item for item in aav2.features if item.kind == "repeat_region"]
        assert len(repeats) == 2
        assert repeats[0].start == 0
        assert repeats[-1].end == 4679
        assert all(item.length == 145 for item in repeats)

    def test_the_spliced_coding_sequences_keep_both_exons(self, aav2) -> None:
        spliced = [item for item in aav2.features if len(item.parts) > 1]
        assert len(spliced) == 2
        for item in spliced:
            assert item.parts[-1].end == 2252
            assert item.length == sum(part.length for part in item.parts)

    def test_rep78_is_the_length_of_its_two_exons(self, aav2) -> None:
        rep = next(
            item
            for item in aav2.features
            if len(item.parts) > 1 and item.start == 320
        )
        assert rep.length == 1586 + 25

    def test_nothing_crosses_the_origin_of_a_linear_genome(self, aav2) -> None:
        assert aav2.crossing() == ()

    def test_a_linear_track_draws_every_feature_once(self, aav2) -> None:
        selection = select(aav2, kinds=["CDS", "repeat_region"])
        drawing = draw_linear(aav2, selection, Style(width=900))
        assert drawing.split == ()
        assert drawing.label_collisions == []


class TestCassette:
    def test_the_table_and_a_length_draw_a_cassette(self, examples) -> None:
        from plasmid_map.reader import annotate, blank

        construct, added = annotate(
            blank(2744, "linear", "cassette").first,
            table=str(examples / "cassette.tsv"),
        )
        assert added == 9
        drawing = draw_linear(construct, select(construct), Style(width=820))
        assert drawing.drawn == 9
        assert drawing.label_collisions == []

    def test_the_primer_reads_backwards(self, examples) -> None:
        from plasmid_map.reader import annotate, blank

        construct, _added = annotate(
            blank(2744, "linear").first, table=str(examples / "cassette.tsv")
        )
        primer = next(
            item for item in construct.features if item.kind == "primer_bind"
        )
        assert primer.strand == -1
        assert (primer.start + 1, primer.end) == (1760, 1780)
