#!/usr/bin/env python3
"""
    File name: flatten_annotation.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Flatten a GTF annotation into DEXSeq-style disjoint exonic_part
                 bins with unique exon_id attributes, handling opposite-strand
                 overlaps and multi-locus / multi-chromosome duplicate gene IDs.
"""

import argparse
import collections
import gzip
import os
import sys


DEFAULT_LOCUS_GAP = 250_000
STRAND_RANK = {"+": 0, "-": 1, ".": 2}


class UnionFind:
    """Disjoint-set structure with path compression."""

    def __init__(self, size):
        self.parent = list(range(size))

    def find(self, item):
        parent = self.parent
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(self, left, right):
        root_l = self.find(left)
        root_r = self.find(right)
        if root_l != root_r:
            if root_l < root_r:
                self.parent[root_r] = root_l
            else:
                self.parent[root_l] = root_r


class LocusInstance:
    """Single genomic locus for a sanitized gene_id on one (chrom, strand)."""

    __slots__ = (
        "chrom",
        "strand",
        "raw_gene_id",
        "start_0",
        "end_0",
        "transcripts",
        "locus_id",
    )

    def __init__(self, chrom, strand, raw_gene_id, start_0, end_0, transcripts):
        self.chrom = chrom
        self.strand = strand
        self.raw_gene_id = raw_gene_id
        self.start_0 = start_0
        self.end_0 = end_0
        self.transcripts = transcripts
        self.locus_id = raw_gene_id


def _sanitize_id(raw_value):
    """Normalize characters that break DEXSeq, featureCounts CORE, or GTF quoting."""
    cleaned = raw_value.strip()
    for bad_char in (":", ",", ";", '"', "\t", "\r", "\n", " "):
        if bad_char in cleaned:
            cleaned = cleaned.replace(bad_char, "_")
    return cleaned


def _parse_attr_clause(clause):
    """Parse a single `;`-delimited GTF/GFF attribute clause into `(key, value)`."""
    eq_pos = clause.find("=")
    sp_pos = clause.find(" ")
    tb_pos = clause.find("\t")
    sep_candidates = [pos for pos in (eq_pos, sp_pos, tb_pos) if pos > 0]
    if not sep_candidates:
        return clause, ""

    sep_pos = min(sep_candidates)
    key = clause[:sep_pos].strip()
    rest = clause[sep_pos:].lstrip(" \t=").strip()
    if len(rest) >= 2 and rest[0] == '"' and rest[-1] == '"':
        rest = rest[1:-1].replace('\\"', '"').replace("\\\\", "\\")
    return key, rest.strip()


def _extract_gtf_ids(attr_field):
    """
    Extract `(gene_id, transcript_id)` from a GTF/GFF attribute column,
    ignoring semicolons inside double-quoted attribute values.
    """
    gene_id = ""
    transcript_id = ""
    in_quotes = False
    escape = False
    start = 0
    length = len(attr_field)

    for idx, ch in enumerate(attr_field):
        if escape:
            escape = False
            continue
        if ch == "\\" and in_quotes:
            escape = True
            continue
        if ch == '"':
            in_quotes = not in_quotes
        elif ch == ";" and not in_quotes:
            clause = attr_field[start:idx].strip()
            start = idx + 1
            if clause:
                key, val = _parse_attr_clause(clause)
                if key == "gene_id" and not gene_id:
                    gene_id = val
                elif key == "transcript_id" and not transcript_id:
                    transcript_id = val
                if gene_id and transcript_id:
                    return gene_id, transcript_id

    if start < length:
        clause = attr_field[start:].strip()
        if clause:
            key, val = _parse_attr_clause(clause)
            if key == "gene_id" and not gene_id:
                gene_id = val
            elif key == "transcript_id" and not transcript_id:
                transcript_id = val

    return gene_id, transcript_id


def _merge_intervals(intervals, gap=0):
    """
    Merge 0-based half-open intervals `[start, end)`.
    When `gap is None`, merge strictly overlapping intervals (`start_0 < prev_e`)
    without coalescing book-ended `[a, b)` and `[b, c)` records.
    """
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for start_0, end_0 in ordered[1:]:
        prev_s, prev_e = merged[-1]
        connects = (
            (start_0 < prev_e) if gap is None else (start_0 <= prev_e + gap)
        )
        if connects:
            if end_0 > prev_e:
                merged[-1] = (prev_s, end_0)
        else:
            merged.append((start_0, end_0))
    return merged


def _spans_connect(next_start, active_end, locus_gap):
    """Return whether `next_start` connects to `active_end` given `locus_gap`."""
    if locus_gap > 0:
        return next_start <= active_end + locus_gap
    return next_start < active_end


def _split_transcript_instances(tx_id, exons, raw_tx_spans):
    """
    Build one or more `(tx_id, merged_exons, start_0, end_0)` tuples for a
    `(chrom, strand, gene_id, tx_id)` group, splitting only when multiple
    strictly disjoint `transcript` feature records exist for the same `tx_id`.
    """
    merged_tx_spans = _merge_intervals(raw_tx_spans, gap=None)
    if len(merged_tx_spans) <= 1:
        merged_exons = _merge_intervals(exons, gap=0)
        return (
            [(tx_id, merged_exons, merged_exons[0][0], merged_exons[-1][1])],
            False,
        )

    buckets = [[] for _ in merged_tx_spans]
    for ex_start, ex_end in exons:
        best_idx = None
        best_overlap = 0
        for idx, (ts_0, te_0) in enumerate(merged_tx_spans):
            overlap = min(ex_end, te_0) - max(ex_start, ts_0)
            if overlap > best_overlap:
                best_overlap = overlap
                best_idx = idx
        if best_idx is None:
            best_idx = min(
                range(len(merged_tx_spans)),
                key=lambda idx: max(
                    0,
                    merged_tx_spans[idx][0] - ex_end,
                    ex_start - merged_tx_spans[idx][1],
                ),
            )
        buckets[best_idx].append((ex_start, ex_end))

    instances = []
    for bucket in buckets:
        if not bucket:
            continue
        merged_exons = _merge_intervals(bucket, gap=0)
        instances.append(
            (tx_id, merged_exons, merged_exons[0][0], merged_exons[-1][1])
        )
    return instances, len(instances) > 1


def _build_locus_instances(
    exons_by_csg,
    gene_spans_by_csg,
    tx_spans_by_csgt,
    locus_gap=DEFAULT_LOCUS_GAP,
):
    """
    Partition each `(chrom, strand, gene_id)` into one or more `LocusInstance`s.
    Opposite-strand and cross-chromosome records are always distinct; same-strand
    records split strictly (`gap=0`) when explicit `gene` records or duplicated
    `transcript` records define separate loci, or across gaps exceeding
    `locus_gap` when no `gene` records exist.
    """
    all_instances = []

    for (chrom, strand, gene_id), tx_exons_map in exons_by_csg.items():
        tx_instances = []
        has_dup_tx = False
        for tx_id, raw_exons in tx_exons_map.items():
            raw_tx_spans = tx_spans_by_csgt.get((chrom, strand, gene_id, tx_id), ())
            split_insts, was_split = _split_transcript_instances(
                tx_id, raw_exons, raw_tx_spans
            )
            tx_instances.extend(split_insts)
            if was_split:
                has_dup_tx = True

        if len(tx_instances) == 1:
            _, _, start_0, end_0 = tx_instances[0]
            all_instances.append(
                LocusInstance(chrom, strand, gene_id, start_0, end_0, tx_instances)
            )
            continue

        raw_gene_spans = gene_spans_by_csg.get((chrom, strand, gene_id), ())
        if raw_gene_spans:
            bridging_genes = _merge_intervals(raw_gene_spans, gap=None)
            effective_gap = 0
        else:
            bridging_genes = []
            effective_gap = 0 if has_dup_tx else locus_gap

        n_tx = len(tx_instances)
        combined = [
            (tx_inst[2], tx_inst[3], idx) for idx, tx_inst in enumerate(tx_instances)
        ]
        for g_idx, (gs_0, ge_0) in enumerate(bridging_genes):
            combined.append((gs_0, ge_0, n_tx + g_idx))

        combined.sort(key=lambda item: (item[0], item[1], item[2]))
        uf = UnionFind(n_tx + len(bridging_genes))
        max_end = combined[0][1]
        max_node = combined[0][2]
        for start_0, end_0, node in combined[1:]:
            if _spans_connect(start_0, max_end, effective_gap):
                uf.union(max_node, node)
            if end_0 > max_end:
                max_end = end_0
                max_node = node

        comp_map = collections.defaultdict(list)
        for idx, tx_inst in enumerate(tx_instances):
            comp_map[uf.find(idx)].append(tx_inst)

        for group in comp_map.values():
            loc_start = min(item[2] for item in group)
            loc_end = max(item[3] for item in group)
            all_instances.append(
                LocusInstance(chrom, strand, gene_id, loc_start, loc_end, group)
            )

    return all_instances


def _disambiguate_loci(all_instances, reserved_gene_ids):
    """
    Keep single-locus `gene_id`s unchanged and disambiguate multi-locus `gene_id`s
    deterministically as `<gene_id>_loc1`, `<gene_id>_loc2`, ... without collisions.
    """
    used_ids = set(reserved_gene_ids)
    by_gene = collections.defaultdict(list)
    for inst in all_instances:
        by_gene[inst.raw_gene_id].append(inst)

    disambiguated_genes = 0
    disambiguated_loci = 0
    for raw_gene_id in sorted(by_gene):
        group = by_gene[raw_gene_id]
        if len(group) == 1:
            group[0].locus_id = raw_gene_id
            used_ids.add(raw_gene_id)
            continue

        disambiguated_genes += 1
        disambiguated_loci += len(group)
        group.sort(
            key=lambda inst: (
                inst.chrom,
                STRAND_RANK.get(inst.strand, 99),
                inst.start_0,
                inst.end_0,
            )
        )
        loc_num = 1
        for inst in group:
            while True:
                candidate = f"{raw_gene_id}_loc{loc_num}"
                loc_num += 1
                if candidate not in used_ids:
                    break
            used_ids.add(candidate)
            inst.locus_id = candidate

    return disambiguated_genes, disambiguated_loci


def _disjoin_island(island_exons):
    """
    Disjoin a contiguous overlapping island of `(start_0, end_0, locus_id, tx_id)`
    tuples into maximal steps `(seg_start, seg_end, frozenset_of_pairs)`.
    """
    starts_at = collections.defaultdict(list)
    ends_at = collections.defaultdict(list)
    coord_set = set()

    for start_0, end_0, locus_id, tx_id in island_exons:
        pair = (locus_id, tx_id)
        starts_at[start_0].append(pair)
        ends_at[end_0].append(pair)
        coord_set.add(start_0)
        coord_set.add(end_0)

    coords = sorted(coord_set)
    active_counts = collections.Counter()
    steps = []

    for idx in range(len(coords) - 1):
        pos = coords[idx]
        next_pos = coords[idx + 1]

        for pair in ends_at.get(pos, ()):
            remaining = active_counts[pair] - 1
            if remaining <= 0:
                del active_counts[pair]
            else:
                active_counts[pair] = remaining

        for pair in starts_at.get(pos, ()):
            active_counts[pair] += 1

        if active_counts and next_pos > pos:
            active_pairs = frozenset(active_counts.keys())
            if steps and steps[-1][1] == pos and steps[-1][2] == active_pairs:
                steps[-1] = (steps[-1][0], next_pos, active_pairs)
            else:
                steps.append((pos, next_pos, active_pairs))

    return steps


def _flatten_chrom_strand(instances, aggregate_genes=True):
    """
    Flatten all `LocusInstance`s on a single `(chrom, strand)` into a list of
    `(aggregate_id, aggr_start_0, aggr_end_0, parts)` tuples, where each part is
    `(part_start_0, part_end_0, tx_string)`.
    """
    exons = []
    locus_ids = []
    for inst in instances:
        locus_ids.append(inst.locus_id)
        for tx_id, merged_exons, _, _ in inst.transcripts:
            for start_0, end_0 in merged_exons:
                exons.append((start_0, end_0, inst.locus_id, tx_id))

    if not exons:
        return []

    exons.sort(key=lambda item: (item[0], item[1], item[2], item[3]))

    # Partition into contiguous exonic islands on (chrom, strand)
    islands = []
    cur_island = [exons[0]]
    island_end = exons[0][1]
    for exon in exons[1:]:
        if exon[0] < island_end:
            cur_island.append(exon)
            if exon[1] > island_end:
                island_end = exon[1]
        else:
            islands.append(cur_island)
            cur_island = [exon]
            island_end = exon[1]
    islands.append(cur_island)

    if aggregate_genes:
        unique_loci = sorted(set(locus_ids))
        locus_to_idx = {lid: idx for idx, lid in enumerate(unique_loci)}
        uf = UnionFind(len(unique_loci))

        for island in islands:
            first_idx = locus_to_idx[island[0][2]]
            for _, _, lid, _ in island[1:]:
                uf.union(first_idx, locus_to_idx[lid])

        comp_members = collections.defaultdict(list)
        for idx, lid in enumerate(unique_loci):
            comp_members[uf.find(idx)].append(lid)

        locus_to_agg = {}
        for members in comp_members.values():
            agg_id = "+".join(sorted(members))
            for lid in members:
                locus_to_agg[lid] = agg_id

    aggregates = collections.defaultdict(list)
    for island in islands:
        for seg_start, seg_end, active_pairs in _disjoin_island(island):
            if aggregate_genes:
                first_lid = next(iter(active_pairs))[0]
                agg_id = locus_to_agg[first_lid]
            else:
                step_loci = {lid for lid, _ in active_pairs}
                if len(step_loci) > 1:
                    continue
                agg_id = next(iter(step_loci))

            tx_str = "+".join(sorted({tx_id for _, tx_id in active_pairs}))
            aggregates[agg_id].append((seg_start, seg_end, tx_str))

    results = []
    for agg_id, parts in aggregates.items():
        results.append((agg_id, parts[0][0], parts[-1][1], parts))
    return results


def flatten_annotation(
    gtf_path,
    output_path,
    locus_gap=DEFAULT_LOCUS_GAP,
    aggregate_genes=True,
    keep_unstranded=False,
):
    """
    Read `gtf_path` (`.gtf` or `.gtf.gz`) and write a DEXSeq-style flattened GTF
    containing `aggregate_gene` and `exonic_part` features (with `exon_id`) to
    `output_path`.
    """
    if locus_gap < 0:
        raise ValueError(f"locus_gap must be >= 0 (got {locus_gap}).")

    if os.path.abspath(gtf_path) == os.path.abspath(output_path):
        raise ValueError("Input and output GTF paths must differ.")

    allowed_strands = ("+", "-", ".") if keep_unstranded else ("+", "-")
    opener = (
        (lambda path: gzip.open(path, "rt", encoding="utf-8"))
        if str(gtf_path).endswith(".gz")
        else (lambda path: open(path, "r", encoding="utf-8"))
    )

    exons_by_csg = collections.defaultdict(lambda: collections.defaultdict(list))
    gene_spans_by_csg = collections.defaultdict(list)
    tx_spans_by_csgt = collections.defaultdict(list)
    reserved_gene_ids = set()

    skipped_unstranded = 0
    skipped_missing_gene = 0
    skipped_malformed = 0
    fallback_tx_records = 0

    with opener(gtf_path) as handle:
        for raw_line in handle:
            if not raw_line or raw_line[0] == "#":
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                continue

            feat_type = fields[2]
            if feat_type not in ("exon", "gene", "transcript"):
                continue

            chrom = fields[0]
            strand = fields[6]
            try:
                start_1 = int(fields[3])
                end_1 = int(fields[4])
            except ValueError:
                if feat_type == "exon":
                    skipped_malformed += 1
                continue

            start_0 = start_1 - 1
            end_0 = end_1
            if start_1 < 1 or end_0 <= start_0:
                if feat_type == "exon":
                    skipped_malformed += 1
                continue

            raw_gene, raw_tx = _extract_gtf_ids(fields[8])
            gene_id = _sanitize_id(raw_gene) if raw_gene else ""
            if gene_id:
                reserved_gene_ids.add(gene_id)

            if strand not in allowed_strands:
                if feat_type == "exon":
                    skipped_unstranded += 1
                continue

            if not gene_id:
                if feat_type == "exon":
                    skipped_missing_gene += 1
                continue

            if feat_type == "gene":
                gene_spans_by_csg[(chrom, strand, gene_id)].append((start_0, end_0))
                continue

            tx_id = _sanitize_id(raw_tx) if raw_tx else ""
            if feat_type == "transcript":
                if tx_id:
                    tx_spans_by_csgt[(chrom, strand, gene_id, tx_id)].append(
                        (start_0, end_0)
                    )
                continue

            if not tx_id:
                tx_id = gene_id
                fallback_tx_records += 1

            exons_by_csg[(chrom, strand, gene_id)][tx_id].append((start_0, end_0))

    if skipped_unstranded:
        sys.stderr.write(
            f"Skipping {skipped_unstranded} exon record(s) without a usable strand.\n"
        )
    if skipped_missing_gene:
        sys.stderr.write(
            f"Skipping {skipped_missing_gene} exon record(s) missing gene_id.\n"
        )
    if skipped_malformed:
        sys.stderr.write(
            f"Skipping {skipped_malformed} exon record(s) with invalid coordinates.\n"
        )
    if fallback_tx_records:
        sys.stderr.write(
            f"Assigned gene_id as transcript_id for {fallback_tx_records} "
            "exon record(s) missing transcript_id.\n"
        )

    if not exons_by_csg:
        raise ValueError(
            f"No usable 'exon' records with gene_id and allowed strand "
            f"{allowed_strands} found in {gtf_path}."
        )

    all_instances = _build_locus_instances(
        exons_by_csg,
        gene_spans_by_csg,
        tx_spans_by_csgt,
        locus_gap=locus_gap,
    )
    disambig_genes, disambig_loci = _disambiguate_loci(
        all_instances, reserved_gene_ids
    )
    if disambig_genes:
        sys.stderr.write(
            f"Disambiguated {disambig_genes} multi-locus/multi-strand gene_id(s) "
            f"into {disambig_loci} distinct locus instances.\n"
        )

    by_cs = collections.defaultdict(list)
    for inst in all_instances:
        by_cs[(inst.chrom, inst.strand)].append(inst)

    all_aggregates = []
    for (chrom, strand), cs_instances in by_cs.items():
        for agg_id, aggr_start_0, aggr_end_0, parts in _flatten_chrom_strand(
            cs_instances, aggregate_genes=aggregate_genes
        ):
            all_aggregates.append(
                (chrom, aggr_start_0, aggr_end_0, strand, agg_id, parts)
            )

    if not all_aggregates:
        raise ValueError(
            f"No exonic_part features remained after flattening {gtf_path}."
        )

    all_aggregates.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
            STRAND_RANK.get(item[3], 99),
            item[4],
        )
    )

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    source_name = "flatten_annotation.py"
    used_exon_ids = set()

    with open(output_path, "w", encoding="utf-8", newline="\n") as out_fh:
        for chrom, aggr_start_0, aggr_end_0, strand, agg_id, parts in all_aggregates:
            out_fh.write(
                f"{chrom}\t{source_name}\taggregate_gene\t"
                f"{aggr_start_0 + 1}\t{aggr_end_0}\t.\t{strand}\t.\t"
                f'gene_id "{agg_id}";\n'
            )
            num_width = max(3, len(str(len(parts))))
            for part_idx, (part_start_0, part_end_0, tx_str) in enumerate(parts, 1):
                part_num = f"{part_idx:0{num_width}d}"
                base_exon_id = f"E{agg_id}{part_num}"
                exon_id = base_exon_id
                if exon_id in used_exon_ids:
                    dup_idx = 2
                    while f"{base_exon_id}_{dup_idx}" in used_exon_ids:
                        dup_idx += 1
                    exon_id = f"{base_exon_id}_{dup_idx}"
                used_exon_ids.add(exon_id)

                out_fh.write(
                    f"{chrom}\t{source_name}\texonic_part\t"
                    f"{part_start_0 + 1}\t{part_end_0}\t.\t{strand}\t.\t"
                    f'gene_id "{agg_id}"; transcripts "{tx_str}"; '
                    f'exonic_part_number "{part_num}"; exon_id "{exon_id}";\n'
                )


def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Flatten a GTF annotation into DEXSeq-style disjoint exonic_part "
            "bins with unique exon_id attributes."
        )
    )
    parser.add_argument(
        "-g",
        "--gtf",
        required=True,
        help="Input GTF annotation path (.gtf or .gtf.gz).",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output flattened GTF path.",
    )
    parser.add_argument(
        "--locus-gap",
        type=int,
        default=DEFAULT_LOCUS_GAP,
        help=(
            "Maximum gap (bp) when clustering same-strand transcripts sharing a "
            "gene_id without an enclosing gene record "
            f"(default: {DEFAULT_LOCUS_GAP})."
        ),
    )
    parser.add_argument(
        "--no-aggregate-genes",
        action="store_true",
        default=False,
        help=(
            "Do not merge distinct overlapping genes on the same strand; "
            "drop exonic segments shared by multiple genes (DEXSeq -r no style)."
        ),
    )
    parser.add_argument(
        "--keep-unstranded",
        action="store_true",
        default=False,
        help="Retain exon records with strand '.' in a separate unstranded partition.",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    cli_args = parse_args()
    flatten_annotation(
        cli_args.gtf,
        cli_args.output,
        locus_gap=cli_args.locus_gap,
        aggregate_genes=not cli_args.no_aggregate_genes,
        keep_unstranded=cli_args.keep_unstranded,
    )
