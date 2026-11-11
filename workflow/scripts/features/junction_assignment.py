#!/usr/bin/env python3
"""
    File name: junction_assignment.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Assign queryname-grouped BAM alignments to deduplicated splice junctions via STAR jM/jI tags.
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
        "-b", "--bam", type=str, required=True, metavar="in_file.bam",
        help="Queryname-sorted BAM file to process",
    )
    parser.add_argument(
        "-o", "--output", type=str, default=None,
        help="Output CSV or CSV.GZ path (default: <bam_prefix>_junctions.csv)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of htslib decompression threads (default: 1)",
    )
    parser.add_argument(
        "--filterBySJout", dest="filterBySJout", action="store_true", default=True,
        help="Only output junctions included in STAR's SJ.out.tab (jM >= 20; default: True)",
    )
    parser.add_argument(
        "--noFilterBySJout", dest="filterBySJout", action="store_false",
        help="Allow junctions not present in STAR's SJ.out.tab",
    )
    parser.add_argument(
        "--filterNonCanon", action="store_true", default=False,
        help="Exclude non-canonical junctions (jM == 0 or jM == 20)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_name = args.bam[:-4] if args.bam.endswith(".bam") else args.bam.split(".bam")[0]
    output_path = args.output if args.output else (input_name + "_junctions.csv")
    is_gz = output_path.endswith(".gz")

    sj_filter = args.filterBySJout
    nc_filter = args.filterNonCanon

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    opener = (
        (lambda p: gzip.open(p, "wt", newline="", compresslevel=3))
        if is_gz
        else (lambda p: open(p, "w", newline=""))
    )

    header = ["qname", "junction_start", "junction_end"]

    with pysam.AlignmentFile(args.bam, "rb", threads=max(1, args.threads)) as samfile, opener(
        output_path
    ) as out_fh:
        wr = csv.writer(out_fh)
        wr.writerow(header)

        for qname, group in itertools.groupby(samfile, key=lambda r: r.query_name):
            junctions = set()

            for read in group:
                if read.is_unmapped:
                    continue

                has_jm = read.has_tag("jM")
                has_ji = read.has_tag("jI")
                if not (has_jm and has_ji):
                    if read.cigarstring and "N" in read.cigarstring:
                        raise RuntimeError(
                            f"Spliced alignment '{qname}' is missing STAR jM/jI tags. "
                            "Ensure STAR was run with '--outSAMattributes ... jI jM' when features.junctions is enabled."
                        )
                    continue

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
                    start = int(ji[2 * k])
                    end = int(ji[2 * k + 1])
                    junctions.add((start, end))

            if junctions:
                ordered = sorted(junctions)
                starts = "+".join(str(s) for s, _ in ordered)
                ends = "+".join(str(e) for _, e in ordered)
                wr.writerow([qname, starts, ends])


if __name__ == "__main__":
    main()
