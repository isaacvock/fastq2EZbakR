#!/usr/bin/env python3
"""
    File name: extract_splice_sites.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Standalone replacement for HISAT2's hisat2_extract_splice_sites.py
                 (omitted from Bioconda's hisat-3n package). Emits unique 0-based
                 splice sites (`chrom\tleft_0\tright_0\tstrand`) from .gtf/.gtf.gz.
"""

import collections
import gzip
import sys


def _extract_transcript_id(attr_field):
    """Extract `transcript_id` from a GTF/GFF attribute column (quote-aware)."""
    in_quotes = False
    escape = False
    start = 0
    length = len(attr_field)

    def _check_clause(clause):
        if not clause:
            return ""
        eq_pos = clause.find("=")
        sp_pos = clause.find(" ")
        tb_pos = clause.find("\t")
        seps = [pos for pos in (eq_pos, sp_pos, tb_pos) if pos > 0]
        if not seps:
            return ""
        sep_pos = min(seps)
        key = clause[:sep_pos].strip()
        if key != "transcript_id":
            return ""
        val = clause[sep_pos:].lstrip(" \t=").strip()
        if len(val) >= 2 and val[0] == '"' and val[-1] == '"':
            val = val[1:-1].replace('\\"', '"').replace("\\\\", "\\")
        return val.strip()

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
            found = _check_clause(attr_field[start:idx].strip())
            if found:
                return found
            start = idx + 1

    if start < length:
        return _check_clause(attr_field[start:].strip())
    return ""


def extract_splice_sites(gtf_path):
    """
    Return sorted list of unique `(chrom, left_0, right_0, strand)` splice sites
    matching `hisat2_extract_splice_sites.py` (merging <=5 bp exon gaps).
    """
    opener = (
        (lambda p: gzip.open(p, "rt", encoding="utf-8"))
        if str(gtf_path).endswith(".gz")
        else (lambda p: open(p, "r", encoding="utf-8"))
    )

    trans_exons = collections.defaultdict(list)
    with opener(gtf_path) as handle:
        for raw_line in handle:
            if not raw_line or raw_line[0] == "#":
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            strand = fields[6]
            if strand not in ("+", "-"):
                continue
            try:
                left_1 = int(fields[3])
                right_1 = int(fields[4])
            except ValueError:
                continue
            if left_1 < 1 or right_1 < left_1:
                continue
            tx_id = _extract_transcript_id(fields[8])
            if not tx_id:
                continue
            trans_exons[(fields[0], strand, tx_id)].append((left_1, right_1))

    junctions = set()
    for (chrom, strand, _), exons in trans_exons.items():
        if len(exons) < 2:
            continue
        exons.sort()
        merged = [exons[0]]
        for left_1, right_1 in exons[1:]:
            if left_1 - merged[-1][1] <= 5:
                if right_1 > merged[-1][1]:
                    merged[-1] = (merged[-1][0], right_1)
            else:
                merged.append((left_1, right_1))

        for idx in range(1, len(merged)):
            junctions.add((chrom, merged[idx - 1][1] - 1, merged[idx][0] - 1, strand))

    return sorted(junctions)


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        sys.stderr.write("Usage: extract_splice_sites.py <annotation.gtf[.gz]>\n")
        sys.exit(1)

    for chrom, left_0, right_0, strand in extract_splice_sites(args[0]):
        sys.stdout.write(f"{chrom}\t{left_0}\t{right_0}\t{strand}\n")


if __name__ == "__main__":
    main()
