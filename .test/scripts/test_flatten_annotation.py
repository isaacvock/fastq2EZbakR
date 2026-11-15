#!/usr/bin/env python3
"""Regression and edge-case test suite for flatten_annotation.py."""

import collections
import gzip
import importlib.util
import io
import os
import random
import subprocess
import sys
import tempfile
from contextlib import redirect_stderr


def _load_flattener(module_path):
    spec = importlib.util.spec_from_file_location("flatten_annotation", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _gtf_row(
    start,
    end,
    strand="+",
    gene="G1",
    tx="T1",
    feat_type="exon",
    chrom="chr1",
    extra="",
):
    attrs = []
    if gene is not None:
        attrs.append(f'gene_id "{gene}";')
    if tx is not None:
        attrs.append(f'transcript_id "{tx}";')
    if extra:
        attrs.append(extra)
    attr_str = " ".join(attrs)
    return f"{chrom}\tfixture\t{feat_type}\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}"


def _parse_flat_gtf(path):
    records = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            assert not line.endswith("\r\n"), "Expected LF line endings"
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 9, f"Expected 9 GTF columns, got {len(fields)}"
            chrom, source, feat_type, start_s, end_s, score, strand, frame, attr_s = (
                fields
            )
            start_1 = int(start_s)
            end_1 = int(end_s)
            assert 1 <= start_1 <= end_1, f"Invalid interval {start_1}..{end_1}"
            assert score == "." and frame == "."

            attrs = {}
            for raw_clause in attr_s.strip().split(";"):
                clause = raw_clause.strip()
                if not clause:
                    continue
                key, val = clause.split(" ", 1)
                assert val.startswith('"') and val.endswith('"'), f"Unquoted attr {val}"
                attrs[key] = val[1:-1]

            records.append(
                {
                    "chrom": chrom,
                    "source": source,
                    "type": feat_type,
                    "start": start_1,
                    "end": end_1,
                    "strand": strand,
                    "attrs": attrs,
                }
            )
    return records


def _assert_structural_invariants(records):
    aggr_ids = set()
    exon_ids = set()
    cur_aggr = None
    prev_part_end = None

    for rec in records:
        if rec["type"] == "aggregate_gene":
            gid = rec["attrs"]["gene_id"]
            assert gid not in aggr_ids, f"Duplicate aggregate_gene id: {gid}"
            aggr_ids.add(gid)
            cur_aggr = rec
            prev_part_end = None
        elif rec["type"] == "exonic_part":
            assert cur_aggr is not None, "exonic_part appeared before aggregate_gene"
            assert rec["chrom"] == cur_aggr["chrom"]
            assert rec["strand"] == cur_aggr["strand"]
            assert rec["attrs"]["gene_id"] == cur_aggr["attrs"]["gene_id"]
            assert cur_aggr["start"] <= rec["start"] <= rec["end"] <= cur_aggr["end"]
            if prev_part_end is not None:
                assert prev_part_end < rec["start"], (
                    f"Non-disjoint or unsorted exonic_part: {prev_part_end} vs "
                    f"{rec['start']}"
                )
            prev_part_end = rec["end"]

            eid = rec["attrs"]["exon_id"]
            assert eid not in exon_ids, f"Duplicate exon_id: {eid}"
            exon_ids.add(eid)
        else:
            raise AssertionError(f"Unexpected feature type: {rec['type']}")


def _reference_clean_dexseq(gtf_path):
    """
    Independent coordinate-step reference matching legacy DEXSeq + exon_ID.sh
    (with deterministic sorted `+`-joining) on clean single-locus GTFs.
    """
    cs_coords = collections.defaultdict(set)
    cs_exons = collections.defaultdict(list)

    with open(gtf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, strand = fields[0], fields[6]
            start_0 = int(fields[3]) - 1
            end_0 = int(fields[4])
            attrs = {}
            for item in fields[8].split(";"):
                item = item.strip()
                if not item:
                    continue
                k, v = item.split(" ", 1)
                attrs[k] = v.strip('"')
            gid = attrs["gene_id"].replace(":", "_")
            tid = attrs["transcript_id"]
            cs_coords[(chrom, strand)].update((start_0, end_0))
            cs_exons[(chrom, strand)].append((start_0, end_0, gid, tid))

    steps_by_cs = {}
    gene_sets = collections.defaultdict(set)

    for cs_key, coord_set in cs_coords.items():
        coords = sorted(coord_set)
        idx_map = {pos: i for i, pos in enumerate(coords)}
        step_sets = [set() for _ in range(len(coords) - 1)]
        for s_0, e_0, gid, tid in cs_exons[cs_key]:
            for step_idx in range(idx_map[s_0], idx_map[e_0]):
                step_sets[step_idx].add((gid, tid))

        coalesced = []
        for step_idx, pairs in enumerate(step_sets):
            if not pairs:
                continue
            s_0, e_0 = coords[step_idx], coords[step_idx + 1]
            frozen = frozenset(pairs)
            if coalesced and coalesced[-1][1] == s_0 and coalesced[-1][2] == frozen:
                coalesced[-1] = (coalesced[-1][0], e_0, frozen)
            else:
                coalesced.append((s_0, e_0, frozen))
        steps_by_cs[cs_key] = coalesced

        for _, _, pairs in coalesced:
            full_set = set()
            for gid, _ in pairs:
                full_set.add(gid)
                full_set |= gene_sets[gid]
            for gid in full_set:
                gene_sets[gid] = full_set

    aggregates = collections.defaultdict(list)
    aggr_meta = {}
    for (chrom, strand), coalesced in steps_by_cs.items():
        for s_0, e_0, pairs in coalesced:
            any_gid = next(iter(pairs))[0]
            agg_id = "+".join(sorted(gene_sets[any_gid]))
            tx_str = "+".join(sorted({tid for _, tid in pairs}))
            aggregates[agg_id].append((s_0, e_0, tx_str))
            aggr_meta[agg_id] = (chrom, strand)

    expected = []
    ordered_aggs = sorted(
        aggregates.keys(),
        key=lambda aid: (
            aggr_meta[aid][0],
            aggregates[aid][0][0],
            aggregates[aid][-1][1],
            0 if aggr_meta[aid][1] == "+" else 1,
            aid,
        ),
    )
    for aid in ordered_aggs:
        chrom, strand = aggr_meta[aid]
        parts = aggregates[aid]
        expected.append(
            (
                "aggregate_gene",
                chrom,
                parts[0][0] + 1,
                parts[-1][1],
                strand,
                aid,
                "",
                "",
            )
        )
        width = max(3, len(str(len(parts))))
        for idx, (s_0, e_0, tx_str) in enumerate(parts, 1):
            pnum = f"{idx:0{width}d}"
            expected.append(
                (
                    "exonic_part",
                    chrom,
                    s_0 + 1,
                    e_0,
                    strand,
                    aid,
                    tx_str,
                    f"E{aid}{pnum}",
                )
            )
    return expected


def run_all_tests(flattener_path, test_genome_gtf):
    mod = _load_flattener(flattener_path)

    with tempfile.TemporaryDirectory(prefix="flatten_test_") as tmp_dir:
        # 1. Exact parity & shuffle determinism on .test/data/annotation/genome.gtf
        flat_genome = os.path.join(tmp_dir, "genome_flat.gtf")
        mod.flatten_annotation(test_genome_gtf, flat_genome)
        parsed_genome = _parse_flat_gtf(flat_genome)
        _assert_structural_invariants(parsed_genome)

        observed_tuples = [
            (
                r["type"],
                r["chrom"],
                r["start"],
                r["end"],
                r["strand"],
                r["attrs"]["gene_id"],
                r["attrs"].get("transcripts", ""),
                r["attrs"].get("exon_id", ""),
            )
            for r in parsed_genome
        ]
        expected_tuples = _reference_clean_dexseq(test_genome_gtf)
        assert observed_tuples == expected_tuples, (
            "Mismatch against DEXSeq reference on .test/data/annotation/genome.gtf"
        )

        with open(test_genome_gtf, "r", encoding="utf-8") as fh:
            body_lines = [ln for ln in fh if not ln.startswith("#")]
        rng = random.Random(42)
        rng.shuffle(body_lines)
        shuffled_in = os.path.join(tmp_dir, "shuffled_in.gtf")
        shuffled_out = os.path.join(tmp_dir, "shuffled_out.gtf")
        with open(shuffled_in, "w", encoding="utf-8", newline="\n") as fh:
            fh.writelines(body_lines)
        mod.flatten_annotation(shuffled_in, shuffled_out)
        with open(flat_genome, "r", encoding="utf-8") as fh1, open(
            shuffled_out, "r", encoding="utf-8"
        ) as fh2:
            assert fh1.read() == fh2.read(), "Output depended on input GTF row order"

        # 2a. Opposite-strand overlap with DIFFERENT gene_ids (no cross-strand split)
        in_2a = os.path.join(tmp_dir, "opp_diff.gtf")
        out_2a = os.path.join(tmp_dir, "opp_diff_flat.gtf")
        with open(in_2a, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(
                "\n".join(
                    [
                        _gtf_row(100, 300, "+", "GenePlus", "txP"),
                        _gtf_row(200, 400, "-", "GeneMinus", "txM"),
                    ]
                )
                + "\n"
            )
        mod.flatten_annotation(in_2a, out_2a)
        recs_2a = _parse_flat_gtf(out_2a)
        _assert_structural_invariants(recs_2a)
        parts_2a = [r for r in recs_2a if r["type"] == "exonic_part"]
        assert len(parts_2a) == 2
        assert (
            parts_2a[0]["start"],
            parts_2a[0]["end"],
            parts_2a[0]["strand"],
            parts_2a[0]["attrs"]["exon_id"],
        ) == (100, 300, "+", "EGenePlus001")
        assert (
            parts_2a[1]["start"],
            parts_2a[1]["end"],
            parts_2a[1]["strand"],
            parts_2a[1]["attrs"]["exon_id"],
        ) == (200, 400, "-", "EGeneMinus001")

        # 2b. Opposite-strand overlap with SAME gene_id (StringTie MSTRG.* super-locus)
        #     plus strand-specific overlapping neighbor genes (must NOT bridge strands)
        in_2b = os.path.join(tmp_dir, "opp_same.gtf")
        out_2b = os.path.join(tmp_dir, "opp_same_flat.gtf")
        with open(in_2b, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(
                "\n".join(
                    [
                        _gtf_row(100, 300, "+", "MSTRG.42", "MSTRG.42.1"),
                        _gtf_row(500, 700, "+", "MSTRG.42", "MSTRG.42.1"),
                        _gtf_row(200, 600, "-", "MSTRG.42", "MSTRG.42.2"),
                        _gtf_row(650, 800, "+", "NeighborPlus", "txNP"),
                        _gtf_row(150, 250, "-", "NeighborMinus", "txNM"),
                    ]
                )
                + "\n"
            )
        stderr_buf = io.StringIO()
        with redirect_stderr(stderr_buf):
            mod.flatten_annotation(in_2b, out_2b)
        assert "Disambiguated 1" in stderr_buf.getvalue()
        recs_2b = _parse_flat_gtf(out_2b)
        _assert_structural_invariants(recs_2b)
        aggs_2b = {
            r["attrs"]["gene_id"]: (r["strand"], r["start"], r["end"])
            for r in recs_2b
            if r["type"] == "aggregate_gene"
        }
        assert aggs_2b == {
            "MSTRG.42_loc1+NeighborPlus": ("+", 100, 800),
            "MSTRG.42_loc2+NeighborMinus": ("-", 150, 600),
        }

        # 3a. Same gene_id on different chromosomes (PAR gene; chrY start < chrX end)
        # 3b. Same gene_id on opposite strands (non-overlapping)
        # 3f. Pre-existing `<gene>_loc1` collision avoidance
        in_3ab = os.path.join(tmp_dir, "multi_chr.gtf")
        out_3ab = os.path.join(tmp_dir, "multi_chr_flat.gtf")
        with open(in_3ab, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(
                "\n".join(
                    [
                        _gtf_row(5000, 6000, "+", "PLCXD1", "txX", chrom="chrX"),
                        _gtf_row(1000, 2000, "+", "PLCXD1", "txY", chrom="chrY"),
                        _gtf_row(100, 200, "+", "INV", "txI1", chrom="chr1"),
                        _gtf_row(500, 600, "-", "INV", "txI2", chrom="chr1"),
                        _gtf_row(9000, 9100, "+", "FOO_loc1", "txLit", chrom="chr1"),
                        _gtf_row(2000, 2100, "+", "FOO", "txF1", chrom="chr1"),
                        _gtf_row(100, 200, "+", "FOO", "txF2", chrom="chr2"),
                    ]
                )
                + "\n"
            )
        mod.flatten_annotation(in_3ab, out_3ab)
        recs_3ab = _parse_flat_gtf(out_3ab)
        _assert_structural_invariants(recs_3ab)
        aggs_3ab = {
            r["attrs"]["gene_id"]: (r["chrom"], r["strand"], r["start"], r["end"])
            for r in recs_3ab
            if r["type"] == "aggregate_gene"
        }
        assert aggs_3ab["PLCXD1_loc1"] == ("chrX", "+", 5000, 6000)
        assert aggs_3ab["PLCXD1_loc2"] == ("chrY", "+", 1000, 2000)
        assert aggs_3ab["INV_loc1"] == ("chr1", "+", 100, 200)
        assert aggs_3ab["INV_loc2"] == ("chr1", "-", 500, 600)
        assert aggs_3ab["FOO_loc1"] == ("chr1", "+", 9000, 9100)
        assert aggs_3ab["FOO_loc2"] == ("chr1", "+", 2000, 2100)
        assert aggs_3ab["FOO_loc3"] == ("chr2", "+", 100, 200)

        # 3c, 3d, 3e. Disjoint same-strand loci vs single-locus GENCODE fragment ENSTs
        in_3cde = os.path.join(tmp_dir, "same_strand_loci.gtf")
        out_3cde = os.path.join(tmp_dir, "same_strand_loci_flat.gtf")
        with open(in_3cde, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(
                "\n".join(
                    [
                        # 3c: RNU6 duplicated at 1kb and 500kb without gene rows
                        _gtf_row(1000, 1100, "+", "RNU6", "u6_copy1"),
                        _gtf_row(1050, 1200, "+", "LocalA", "txLA"),
                        _gtf_row(500000, 500100, "+", "RNU6", "u6_copy2"),
                        _gtf_row(500050, 500200, "+", "LocalB", "txLB"),
                        # 3d: Duplicated gene & transcript_id with gene/tx rows
                        _gtf_row(10000, 10500, "+", "DUPG", None, feat_type="gene"),
                        _gtf_row(
                            10000, 10500, "+", "DUPG", "NM_1", feat_type="transcript"
                        ),
                        _gtf_row(10000, 10200, "+", "DUPG", "NM_1"),
                        _gtf_row(10300, 10500, "+", "DUPG", "NM_1"),
                        _gtf_row(90000, 90500, "+", "DUPG", None, feat_type="gene"),
                        _gtf_row(
                            90000, 90500, "+", "DUPG", "NM_1", feat_type="transcript"
                        ),
                        _gtf_row(90000, 90200, "+", "DUPG", "NM_1"),
                        _gtf_row(90300, 90500, "+", "DUPG", "NM_1"),
                        # 3d2: Duplicated transcript_id 5kb apart WITHOUT gene rows
                        _gtf_row(
                            700000,
                            700500,
                            "+",
                            "DUPTX",
                            "NM_2",
                            feat_type="transcript",
                        ),
                        _gtf_row(700000, 700500, "+", "DUPTX", "NM_2"),
                        _gtf_row(
                            705000,
                            705500,
                            "+",
                            "DUPTX",
                            "NM_2",
                            feat_type="transcript",
                        ),
                        _gtf_row(705000, 705500, "+", "DUPTX", "NM_2"),
                        # 3e: Single-locus GENCODE gene >250kb enclosing fragment ENSTs
                        _gtf_row(
                            20000, 450000, "+", "ENSG001", None, feat_type="gene"
                        ),
                        _gtf_row(20000, 21000, "+", "ENSG001", "ENST1"),
                        _gtf_row(448000, 450000, "+", "ENSG001", "ENST2"),
                    ]
                )
                + "\n"
            )
        mod.flatten_annotation(in_3cde, out_3cde)
        recs_3cde = _parse_flat_gtf(out_3cde)
        _assert_structural_invariants(recs_3cde)
        aggs_3cde = {
            r["attrs"]["gene_id"]: (r["start"], r["end"])
            for r in recs_3cde
            if r["type"] == "aggregate_gene"
        }
        assert aggs_3cde["LocalA+RNU6_loc1"] == (1000, 1200)
        assert aggs_3cde["LocalB+RNU6_loc2"] == (500000, 500200)
        assert aggs_3cde["DUPG_loc1"] == (10000, 10500)
        assert aggs_3cde["DUPG_loc2"] == (90000, 90500)
        assert aggs_3cde["DUPTX_loc1"] == (700000, 700500)
        assert aggs_3cde["DUPTX_loc2"] == (705000, 705500)
        assert aggs_3cde["ENSG001"] == (20000, 450000)

        # 3g. Custom locus_gap=0 splits short-gap isoforms lacking gene records
        in_3g = os.path.join(tmp_dir, "short_gap.gtf")
        out_3g_def = os.path.join(tmp_dir, "short_gap_default.gtf")
        out_3g_zero = os.path.join(tmp_dir, "short_gap_zero.gtf")
        with open(in_3g, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(
                "\n".join(
                    [
                        _gtf_row(1000, 1500, "+", "SHORTGAP", "txSG1"),
                        _gtf_row(6000, 6500, "+", "SHORTGAP", "txSG2"),
                    ]
                )
                + "\n"
            )
        mod.flatten_annotation(in_3g, out_3g_def)
        aggs_3g_def = [
            r["attrs"]["gene_id"]
            for r in _parse_flat_gtf(out_3g_def)
            if r["type"] == "aggregate_gene"
        ]
        assert aggs_3g_def == ["SHORTGAP"]

        mod.flatten_annotation(in_3g, out_3g_zero, locus_gap=0)
        aggs_3g_zero = [
            r["attrs"]["gene_id"]
            for r in _parse_flat_gtf(out_3g_zero)
            if r["type"] == "aggregate_gene"
        ]
        assert aggs_3g_zero == ["SHORTGAP_loc1", "SHORTGAP_loc2"]

        # 4a, 4b. Unstranded exons, quoted semicolons, colon sanitization,
        #         book-ended exons, missing transcript_id, and .gtf.gz input
        in_4_gz = os.path.join(tmp_dir, "quirky.gtf.gz")
        out_4_skip = os.path.join(tmp_dir, "quirky_skip.gtf")
        out_4_keep = os.path.join(tmp_dir, "quirky_keep.gtf")
        quirky_lines = [
            "##gff-version 2",
            (
                'chr1\tsrc\texon\t101\t200\t.\t+\t.\t'
                'note "has; semicolon; and \\"transcript_id \\"FAKE\\"\\""; '
                'gene_id "HGNC:99"; transcript_id "TX_A"'
            ),
            (
                'chr1\tsrc\texon\t201\t300\t.\t+\t.\t'
                'gene_id "HGNC:99"; transcript_id "TX_A";'
            ),
            (
                'chr1\tsrc\texon\t150\t250\t.\t.\t.\t'
                'gene_id "HGNC:99"; transcript_id "TX_UNSTR";'
            ),
            "chr1\tsrc\texon\t500\t600\t.\t-\t.\tgene_id=NO_TX_GENE",
            (
                'chr1\tsrc\texon\t700\t800\t.\t+\t.\t'
                'note "missing gene id"; transcript_id "ORPHAN";'
            ),
        ]
        with gzip.open(in_4_gz, "wt", encoding="utf-8", newline="\n") as gz_fh:
            gz_fh.write("\n".join(quirky_lines) + "\n")

        buf_4 = io.StringIO()
        with redirect_stderr(buf_4):
            mod.flatten_annotation(in_4_gz, out_4_skip, keep_unstranded=False)
        msg_4 = buf_4.getvalue()
        assert "Skipping 1 exon record(s) without a usable strand" in msg_4
        assert "Skipping 1 exon record(s) missing gene_id" in msg_4
        assert "Assigned gene_id as transcript_id for 1" in msg_4

        recs_4_skip = _parse_flat_gtf(out_4_skip)
        _assert_structural_invariants(recs_4_skip)
        parts_4_skip = [r for r in recs_4_skip if r["type"] == "exonic_part"]
        assert len(parts_4_skip) == 2
        # Notice HGNC_99 is NOT suffixed _loc1 when its '.' isoform was skipped!
        assert parts_4_skip[0]["attrs"]["gene_id"] == "HGNC_99"
        assert (parts_4_skip[0]["start"], parts_4_skip[0]["end"]) == (101, 300)
        assert parts_4_skip[0]["attrs"]["transcripts"] == "TX_A"
        assert parts_4_skip[1]["attrs"]["gene_id"] == "NO_TX_GENE"
        assert parts_4_skip[1]["attrs"]["transcripts"] == "NO_TX_GENE"

        mod.flatten_annotation(in_4_gz, out_4_keep, keep_unstranded=True)
        recs_4_keep = _parse_flat_gtf(out_4_keep)
        _assert_structural_invariants(recs_4_keep)
        aggs_4_keep = {
            r["attrs"]["gene_id"]: r["strand"]
            for r in recs_4_keep
            if r["type"] == "aggregate_gene"
        }
        assert aggs_4_keep["HGNC_99_loc1"] == "+"
        assert aggs_4_keep["HGNC_99_loc2"] == "."

        # 4c. --no-aggregate-genes drops multi-gene overlap while keeping flanks
        in_4c = os.path.join(tmp_dir, "no_agg.gtf")
        out_4c = os.path.join(tmp_dir, "no_agg_flat.gtf")
        with open(in_4c, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(
                "\n".join(
                    [
                        _gtf_row(101, 400, "+", "GeneA", "txA"),
                        _gtf_row(201, 300, "+", "GeneB", "txB"),
                    ]
                )
                + "\n"
            )
        mod.flatten_annotation(in_4c, out_4c, aggregate_genes=False)
        recs_4c = _parse_flat_gtf(out_4c)
        _assert_structural_invariants(recs_4c)
        parts_4c = [r for r in recs_4c if r["type"] == "exonic_part"]
        assert [(p["start"], p["end"], p["attrs"]["exon_id"]) for p in parts_4c] == [
            (101, 200, "EGeneA001"),
            (301, 400, "EGeneA002"),
        ]

        # 4d. >999 parts dynamic zero-padding and exon_id collision guard
        #     (GENE1 part 1001 -> EGENE11001 vs GENE11 part 001 -> EGENE11001_2)
        in_4d = os.path.join(tmp_dir, "many_parts.gtf")
        out_4d = os.path.join(tmp_dir, "many_parts_flat.gtf")
        with open(in_4d, "w", encoding="utf-8", newline="\n") as fh:
            for idx in range(1001):
                s_1 = idx * 10 + 1
                e_1 = s_1 + 4
                fh.write(_gtf_row(s_1, e_1, "+", "GENE1", "txBig") + "\n")
            fh.write(_gtf_row(50000, 50100, "+", "GENE11", "txSmall") + "\n")
        mod.flatten_annotation(in_4d, out_4d)
        recs_4d = _parse_flat_gtf(out_4d)
        _assert_structural_invariants(recs_4d)
        parts_4d = [r for r in recs_4d if r["type"] == "exonic_part"]
        assert parts_4d[0]["attrs"]["exon_id"] == "EGENE10001"
        assert parts_4d[1000]["attrs"]["exon_id"] == "EGENE11001"
        assert parts_4d[1001]["attrs"]["exon_id"] == "EGENE11001_2"

        # 4e. Error guardrails & CLI subprocess execution
        try:
            mod.flatten_annotation(in_2a, in_2a)
            raise AssertionError("Expected ValueError when input == output")
        except ValueError as exc:
            assert "differ" in str(exc)

        empty_in = os.path.join(tmp_dir, "no_exons.gtf")
        with open(empty_in, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(_gtf_row(1, 100, "+", "G1", None, feat_type="gene") + "\n")
        try:
            mod.flatten_annotation(empty_in, os.path.join(tmp_dir, "unused.gtf"))
            raise AssertionError("Expected ValueError when no usable exons exist")
        except ValueError as exc:
            assert "No usable" in str(exc)

        cli_out = os.path.join(tmp_dir, "cli_out.gtf")
        proc = subprocess.run(
            [sys.executable, flattener_path, "-g", in_2b, "-o", cli_out],
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, f"CLI failed: {proc.stderr}"
        assert len(_parse_flat_gtf(cli_out)) == len(recs_2b)


def main():
    if len(sys.argv) != 4:
        sys.stderr.write(
            "Usage: test_flatten_annotation.py <flatten_annotation.py> "
            "<genome.gtf> <sentinel.ok>\n"
        )
        sys.exit(1)

    flattener_path, genome_gtf, sentinel_out = sys.argv[1:4]
    run_all_tests(flattener_path, genome_gtf)

    out_dir = os.path.dirname(os.path.abspath(sentinel_out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(sentinel_out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("ok\n")
    print("All flatten_annotation regression tests passed.")


if __name__ == "__main__":
    main()
