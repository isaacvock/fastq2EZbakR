#!/usr/bin/env python3
"""
    File name: transcript_assignment.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Create table of strand-aware transcript equivalence classes (TECs)
                 from a STAR or RSEM transcriptome BAM file.
"""

import csv
import itertools
import sys

import pysam


def normalize_strandedness(strandedness):
    norm = str(strandedness).strip().lower()
    mapping = {
        "reverse": "reverse",
        "antisense": "reverse",
        "r": "reverse",
        "yes": "forward",
        "forward": "forward",
        "sense": "forward",
        "f": "forward",
        "no": "unstranded",
        "none": "unstranded",
        "unstranded": "unstranded",
    }
    if norm not in mapping:
        raise ValueError(
            f"Unsupported strandedness '{strandedness}'. Expected 'reverse', 'yes'/'forward', or 'no'."
        )
    return mapping[norm]


def is_sense_transcriptome_alignment(read, norm_strand):
    """
    Check whether an alignment against a 5'->3' sense transcriptome reference
    matches the library strandedness protocol.
    """
    if norm_strand == "unstranded":
        return True

    if read.is_paired:
        if read.is_read1:
            return read.is_reverse if norm_strand == "reverse" else (not read.is_reverse)
        if read.is_read2:
            return (not read.is_reverse) if norm_strand == "reverse" else read.is_reverse
        return False

    return read.is_reverse if norm_strand == "reverse" else (not read.is_reverse)


def run_transcript_assignment(input_bam, output_csv, strandedness="no", threads=1):
    norm_strand = normalize_strandedness(strandedness)
    header = ["qname", "bamfile_transcripts"]

    with pysam.AlignmentFile(input_bam, "rb", threads=max(1, threads)) as samfile, open(
        output_csv, "w", newline="\n", encoding="utf-8"
    ) as out_fh:
        wr = csv.writer(out_fh, lineterminator="\n")
        wr.writerow(header)

        for qname, group in itertools.groupby(samfile, key=lambda r: r.query_name):
            transcripts = set()
            incompatible_tx = set()

            for r in group:
                if r.is_unmapped or r.reference_name is None:
                    continue
                if r.has_tag("ZW") and r.get_tag("ZW") <= 0:
                    continue

                tx_name = str(r.reference_name)
                if is_sense_transcriptome_alignment(r, norm_strand):
                    transcripts.add(tx_name)
                else:
                    incompatible_tx.add(tx_name)

            if incompatible_tx:
                transcripts.difference_update(incompatible_tx)

            if transcripts:
                wr.writerow([qname, "+".join(sorted(transcripts))])


if "snakemake" in globals():
    with open(snakemake.log[0], "w", encoding="utf-8") as log_fh:
        sys.stdout = log_fh
        sys.stderr = log_fh
        print("Started analysis")
        run_transcript_assignment(
            snakemake.input.bam,
            snakemake.output.table,
            strandedness=snakemake.params.get("strandedness", "no"),
            threads=getattr(snakemake, "threads", 1),
        )
        print("Ended analysis")
elif __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Assign reads to strand-compatible transcripts from a transcriptome BAM"
    )
    parser.add_argument("-b", "--bam", required=True, help="Transcriptome BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV path")
    parser.add_argument(
        "-s",
        "--strandedness",
        default="no",
        help="Library strandedness: 'reverse', 'yes'/'forward', or 'no' (default: no)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Decompression threads"
    )
    cli_args = parser.parse_args()
    run_transcript_assignment(
        cli_args.bam,
        cli_args.output,
        strandedness=cli_args.strandedness,
        threads=cli_args.threads,
    )
