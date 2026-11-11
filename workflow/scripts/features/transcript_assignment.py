"""
    File name: transcript_assignment.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Create table of transcript equivalence classes (TECs) each read mapped to.
"""

import csv
import itertools
import sys

import pysam


def run_transcript_assignment(input_bam, output_csv, threads=1):
    header = ["qname", "bamfile_transcripts"]
    with pysam.AlignmentFile(input_bam, "rb", threads=max(1, threads)) as samfile, open(
        output_csv, "w", newline=""
    ) as out_fh:
        wr = csv.writer(out_fh)
        wr.writerow(header)

        for qname, group in itertools.groupby(samfile, key=lambda r: r.query_name):
            transcripts = set()
            for r in group:
                if r.is_unmapped or r.reference_name is None:
                    continue
                if r.has_tag("ZW"):
                    if r.get_tag("ZW") > 0:
                        transcripts.add(str(r.reference_name))
                else:
                    transcripts.add(str(r.reference_name))

            if transcripts:
                wr.writerow([qname, "+".join(sorted(transcripts))])


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log_fh:
        sys.stdout = log_fh
        sys.stderr = log_fh
        print("Started analysis")
        run_transcript_assignment(
            snakemake.input.bam,
            snakemake.output.table,
            threads=getattr(snakemake, "threads", 1),
        )
        print("Ended analysis")
elif __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Assign reads to transcripts from a transcriptome BAM")
    parser.add_argument("-b", "--bam", required=True, help="Transcriptome BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV path")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Decompression threads")
    cli_args = parser.parse_args()
    run_transcript_assignment(cli_args.bam, cli_args.output, threads=cli_args.threads)
