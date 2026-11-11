#!/usr/bin/env python3
"""
    File name: junction_assignment.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Assign queryname-grouped BAM alignments to deduplicated splice junctions
                 via STAR jM/jI tags when present, with exact CIGAR-N intron fallback
                 for non-STAR aligners (HISAT-3N, HISAT2) and external BAMs.
"""

import argparse
import csv
import gzip
import itertools
import os

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assign reads/read-pairs to splice junctions they overlap"
    )
    required_named = parser.add_argument_group("required named arguments")
    required_named.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        metavar="in_file.bam",
        help="Queryname-sorted BAM file to process",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output CSV or CSV.GZ path (default: <bam_prefix>_junctions.csv)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of htslib decompression threads (default: 1)",
    )
    parser.add_argument(
        "--filterBySJout",
        dest="filterBySJout",
        action="store_true",
        default=True,
        help="When STAR jM tags are present, only output junctions in SJ.out.tab (jM >= 20; default: True)",
    )
    parser.add_argument(
        "--noFilterBySJout",
        dest="filterBySJout",
        action="store_false",
        help="Allow junctions not present in STAR's SJ.out.tab",
    )
    parser.add_argument(
        "--filterNonCanon",
        action="store_true",
        default=False,
        help="When STAR jM tags are present, exclude non-canonical junctions (jM == 0 or jM == 20)",
    )
    parser.add_argument(
        "-g",
        "--gtf",
        type=str,
        default="",
        help="Optional GTF annotation used only if --annotatedOnly is enabled",
    )
    parser.add_argument(
        "--annotatedOnly",
        action="store_true",
        default=False,
        help="Restrict junctions to introns present in --gtf",
    )
    return parser.parse_args()


def _extract_transcript_id(attr_field):
    for raw_item in attr_field.split(";"):
        item = raw_item.strip()
        if not item:
            continue
        if item.startswith("transcript_id"):
            rest = item[len("transcript_id") :].strip()
            if rest.startswith("="):
                rest = rest[1:].strip()
            return rest.strip('"')
    return ""


def load_annotated_introns(gtf_path):
    """Return dict mapping `chrom -> set((intron_start_1based, intron_end_1based))` from GTF exons."""
    opener = (
        (lambda p: gzip.open(p, "rt", encoding="utf-8"))
        if str(gtf_path).endswith(".gz")
        else (lambda p: open(p, "r", encoding="utf-8"))
    )
    exons_by_tx = {}
    with opener(gtf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            tx_id = _extract_transcript_id(fields[8])
            if not tx_id:
                continue
            start_1 = int(fields[3])
            end_1 = int(fields[4])
            if end_1 >= start_1:
                exons_by_tx.setdefault((fields[0], tx_id), []).append((start_1, end_1))

    introns_by_chrom = {}
    for (chrom, _tx_id), ex_list in exons_by_tx.items():
        ex_list.sort()
        merged = [ex_list[0]]
        for s, e in ex_list[1:]:
            ps, pe = merged[-1]
            if s <= pe + 1:
                if e > pe:
                    merged[-1] = (ps, e)
            else:
                merged.append((s, e))
        if len(merged) > 1:
            c_set = introns_by_chrom.setdefault(chrom, set())
            for k in range(len(merged) - 1):
                c_set.add((merged[k][1] + 1, merged[k + 1][0] - 1))
    return introns_by_chrom


def _extract_cigar_introns_1based(read):
    """Extract 1-based closed `(intron_start, intron_end)` pairs from CIGAR `N` operations."""
    cigartuples = read.cigartuples
    if not cigartuples:
        return ()

    r_pos = read.reference_start  # 0-based
    found = []
    for op, length in cigartuples:
        if op == 3:  # N
            found.append((r_pos + 1, r_pos + length))
            r_pos += length
        elif op in (0, 2, 7, 8):  # M, D, =, X
            r_pos += length
    return found


def main():
    args = parse_args()

    input_name = (
        args.bam[:-4] if args.bam.endswith(".bam") else args.bam.split(".bam")[0]
    )
    output_path = args.output if args.output else (input_name + "_junctions.csv")
    is_gz = output_path.endswith(".gz")

    sj_filter = args.filterBySJout
    nc_filter = args.filterNonCanon

    annotated_introns = None
    if args.annotatedOnly:
        if not args.gtf:
            raise ValueError("--annotatedOnly requires -g / --gtf annotation file.")
        annotated_introns = load_annotated_introns(args.gtf)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    opener = (
        (lambda p: gzip.open(p, "wt", newline="\n", compresslevel=3, encoding="utf-8"))
        if is_gz
        else (lambda p: open(p, "w", newline="\n", encoding="utf-8"))
    )

    header = ["qname", "junction_start", "junction_end"]

    with pysam.AlignmentFile(
        args.bam, "rb", threads=max(1, args.threads)
    ) as samfile, opener(output_path) as out_fh:
        wr = csv.writer(out_fh, lineterminator="\n")
        wr.writerow(header)

        for qname, group in itertools.groupby(samfile, key=lambda r: r.query_name):
            junctions = set()

            for read in group:
                if read.is_unmapped:
                    continue

                chrom_annot = (
                    annotated_introns.get(read.reference_name, set())
                    if annotated_introns is not None
                    else None
                )

                if read.has_tag("jM") and read.has_tag("jI"):
                    jm = read.get_tag("jM")
                    ji = read.get_tag("jI")
                    n_junc = len(jm)

                    for k in range(n_junc):
                        motif = jm[k]
                        if motif == -1:
                            continue
                        if sj_filter and motif < 20:
                            continue
                        if nc_filter and (motif == 0 or motif == 20):
                            continue
                        jpair = (int(ji[2 * k]), int(ji[2 * k + 1]))
                        if chrom_annot is not None and jpair not in chrom_annot:
                            continue
                        junctions.add(jpair)
                elif read.cigarstring and "N" in read.cigarstring:
                    for jpair in _extract_cigar_introns_1based(read):
                        if chrom_annot is not None and jpair not in chrom_annot:
                            continue
                        junctions.add(jpair)

            if junctions:
                ordered = sorted(junctions)
                starts = "+".join(str(s) for s, _ in ordered)
                ends = "+".join(str(e) for _, e in ordered)
                wr.writerow([qname, starts, ends])


if __name__ == "__main__":
    main()
