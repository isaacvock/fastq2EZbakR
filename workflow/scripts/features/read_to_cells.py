"""Corrected identities for the same query-sorted alignments used for counting."""

import csv
import pysam


def extract(bam, output):
    previous = None
    missing_cb = missing_ub = total = 0
    with pysam.AlignmentFile(bam, "rb") as alignments, open(output, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["qname", "cell_barcode", "umi"])
        for read in alignments:
            if read.is_paired:
                raise ValueError(f"Single-cell mode requires unpaired cDNA alignments: {read.query_name}")
            if read.query_name == previous:
                raise ValueError(f"Repeated primary single-cell read identity: {read.query_name}")
            previous = read.query_name
            identities = []
            for tag in ("CB", "UB"):
                value = str(read.get_tag(tag)) if read.has_tag(tag) else ""
                identities.append("" if value in ("", "-") else value)
            writer.writerow([read.query_name, *identities])
            total += 1
            missing_cb += not identities[0]
            missing_ub += not identities[1]
    return total, missing_cb, missing_ub


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log:
        total, missing_cb, missing_ub = extract(snakemake.input.samples, snakemake.output.table)
        print(f"Retained {total} reads; missing corrected CB: {missing_cb}; "
              f"missing corrected UB: {missing_ub}. No raw-tag substitution.", file=log)
