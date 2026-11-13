#!/usr/bin/env python3
"""
    File name: call_snps.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Lightweight multi-core mutational-anomaly / SNP caller using samtools mpileup
                 with depth, minimum alt count, allele-frequency, and exact binomial p-value filters.
"""

import argparse
import csv
import datetime
import math
import multiprocessing as mp
import os
import shlex
import subprocess
import sys
import tempfile

import pysam


ALT_BASES = ("A", "C", "G", "T")
VALID_REFS = frozenset(ALT_BASES)


def binom_upper_tail_pval(k, n, p0):
    """
    Compute exact upper-tail binomial p-value P(Binom(n, p0) >= k) in pure Python math.
    Accurate across full float64 range without scipy.
    """
    if k <= 0:
        return 1.0
    if k > n or p0 <= 0.0:
        return 0.0
    if p0 >= 1.0:
        return 1.0

    log_p = math.log(p0)
    log_q = math.log1p(-p0)

    # Standard anomalous site regime: observed rate exceeds null rate
    if k > n * p0:
        log_pmf_k = (
            math.lgamma(n + 1)
            - math.lgamma(k + 1)
            - math.lgamma(n - k + 1)
            + k * log_p
            + (n - k) * log_q
        )
        odds = p0 / (1.0 - p0)
        rel_sum = 1.0
        rel_term = 1.0
        for j in range(k, n):
            rel_term *= ((n - j) / (j + 1)) * odds
            rel_sum += rel_term
            if rel_term <= 1e-15 * rel_sum:
                break
        log_tail = log_pmf_k + math.log(rel_sum)
        if log_tail >= 0.0:
            return 1.0
        return math.exp(log_tail)

    # Fallback when k <= n * p0: sum lower tail j = k-1 down to 0
    j0 = k - 1
    log_pmf_j0 = (
        math.lgamma(n + 1)
        - math.lgamma(j0 + 1)
        - math.lgamma(n - j0 + 1)
        + j0 * log_p
        + (n - j0) * log_q
    )
    inv_odds = (1.0 - p0) / p0
    rel_sum = 1.0
    rel_term = 1.0
    for j in range(j0, 0, -1):
        rel_term *= (j / (n - j + 1)) * inv_odds
        rel_sum += rel_term
        if rel_term <= 1e-15 * rel_sum:
            break
    lower = math.exp(log_pmf_j0 + math.log(rel_sum))
    return max(0.0, min(1.0, 1.0 - lower))


def read_fai_contigs(fai_path):
    """Return list of `(chrom, length)` from a FASTA `.fai` file."""
    contigs = []
    with open(fai_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            length = int(parts[1])
            if length > 0:
                contigs.append((chrom, length))
    return contigs


def get_active_contigs(contigs, bam_paths):
    """
    Filter `.fai` contigs to those with >=1 mapped read in any input BAM using `.bai` metadata.
    Falls back to all `.fai` contigs if index metadata lacks mapped counts despite non-empty BAMs.
    """
    mapped_by_chrom = {}
    any_bam_nonempty = False

    for bam_path in bam_paths:
        with pysam.AlignmentFile(bam_path, "rb") as sf:
            try:
                stats = sf.get_index_statistics()
            except (ValueError, OSError):
                return contigs

            bam_mapped = 0
            for st in stats:
                if st.mapped > 0:
                    mapped_by_chrom[st.contig] = mapped_by_chrom.get(st.contig, 0) + st.mapped
                    bam_mapped += st.mapped

            if bam_mapped > 0:
                any_bam_nonempty = True
            else:
                try:
                    for rec in sf.fetch(until_eof=True):
                        if not rec.is_unmapped:
                            any_bam_nonempty = True
                            break
                except (ValueError, OSError):
                    return contigs

    if not mapped_by_chrom:
        if any_bam_nonempty:
            return contigs
        return []

    return [(c, length) for c, length in contigs if mapped_by_chrom.get(c, 0) > 0]


def build_regions(contigs, threads):
    """Partition active contigs into balanced 1-based closed `(chrom, start, end)` intervals."""
    if not contigs:
        return []

    total_bp = sum(length for _, length in contigs)
    target_chunks = max(1, threads * 4)
    chunk_bp = max(50_000, (total_bp + target_chunks - 1) // target_chunks)

    regions = []
    for chrom, length in contigs:
        start = 1
        while start <= length:
            end = min(length, start + chunk_bp - 1)
            regions.append((chrom, start, end))
            start = end + 1
    return regions


def format_hts_region(chrom, start, end):
    """Format region string for `samtools mpileup -r`, brace-quoting contigs containing colons."""
    if ":" in chrom:
        return f"{{{chrom}}}:{start}-{end}"
    return f"{chrom}:{start}-{end}"


def _scan_region_worker(task):
    (
        region_tuple,
        bam_paths,
        bam_labels,
        fasta_path,
        min_phred,
        min_depth,
        min_alt,
        min_freq,
        err_rate,
        p_cutoff,
        extra_args,
    ) = task

    chrom, start, end = region_tuple
    region_str = format_hts_region(chrom, start, end)
    n_bams = len(bam_paths)

    cmd = [
        "samtools",
        "mpileup",
        "--no-BAQ",
        "--no-output-ins",
        "--no-output-ins",
        "--no-output-del",
        "--no-output-del",
        "--no-output-ends",
        "-Q",
        str(min_phred),
        "-d",
        "100000",
        "-f",
        fasta_path,
        "-r",
        region_str,
        *extra_args,
        *bam_paths,
    ]

    snp_records = []
    expected_min_cols = 3 + 3 * n_bams

    with tempfile.TemporaryFile(mode="w+", encoding="utf-8", errors="replace") as err_fh:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=err_fh,
            text=True,
            encoding="ascii",
            errors="replace",
            bufsize=1024 * 1024,
        )

        assert proc.stdout is not None
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < expected_min_cols:
                continue

            ref = fields[2].upper()
            if ref not in VALID_REFS:
                continue

            # Fast C-level scan per BAM column; skip immediately if zero mismatch bases exist
            has_any_alt = False
            per_bam = []
            pooled_depth = 0
            pooled_A = 0
            pooled_C = 0
            pooled_G = 0
            pooled_T = 0

            for i in range(n_bams):
                raw_dp_str = fields[3 + 3 * i]
                if raw_dp_str == "0":
                    if n_bams > 1:
                        per_bam.append((0, 0, 0, 0, 0))
                    continue

                bases = fields[4 + 3 * i]
                n_ref = bases.count(".") + bases.count(",")
                blen = len(bases)

                if n_ref == blen:
                    pooled_depth += n_ref
                    if n_bams > 1:
                        per_bam.append((n_ref, 0, 0, 0, 0))
                    continue

                n_skip = (
                    bases.count(">")
                    + bases.count("<")
                    + bases.count("*")
                    + bases.count("N")
                    + bases.count("n")
                )
                if n_ref + n_skip == blen:
                    pooled_depth += n_ref
                    if n_bams > 1:
                        per_bam.append((n_ref, 0, 0, 0, 0))
                    continue

                a_cnt = bases.count("A") + bases.count("a")
                c_cnt = bases.count("C") + bases.count("c")
                g_cnt = bases.count("G") + bases.count("g")
                t_cnt = bases.count("T") + bases.count("t")
                dp_i = n_ref + a_cnt + c_cnt + g_cnt + t_cnt

                pooled_depth += dp_i
                if a_cnt or c_cnt or g_cnt or t_cnt:
                    has_any_alt = True
                    pooled_A += a_cnt
                    pooled_C += c_cnt
                    pooled_G += g_cnt
                    pooled_T += t_cnt

                if n_bams > 1:
                    per_bam.append((dp_i, a_cnt, c_cnt, g_cnt, t_cnt))

            if not has_any_alt or pooled_depth < min_depth:
                continue

            pos_int = int(fields[1])
            if pos_int < start or pos_int > end:
                continue
            rec_chrom = fields[0]
            pooled_alts = (("A", pooled_A), ("C", pooled_C), ("G", pooled_G), ("T", pooled_T))

            for alt_idx, (alt_base, p_alt) in enumerate(pooled_alts):
                if alt_base == ref or p_alt < min_alt:
                    continue

                passed = False
                best_dp = pooled_depth
                best_alt = p_alt
                best_freq = p_alt / pooled_depth
                best_pval = 1.0
                best_source = "pooled" if n_bams > 1 else bam_labels[0]

                # 1) Evaluate pooled counts across all input BAMs
                if pooled_depth >= min_depth and p_alt >= min_alt and best_freq >= min_freq:
                    pval_pool = binom_upper_tail_pval(p_alt, pooled_depth, err_rate)
                    if pval_pool <= p_cutoff:
                        passed = True
                        best_pval = pval_pool

                # 2) Evaluate each individual BAM when multiple BAMs are provided
                if not passed and n_bams > 1:
                    alt_col = alt_idx + 1  # (dp_i, A, C, G, T)
                    for bam_idx, stats_i in enumerate(per_bam):
                        dp_i = stats_i[0]
                        if dp_i < min_depth:
                            continue
                        alt_i = stats_i[alt_col]
                        if alt_i < min_alt:
                            continue
                        freq_i = alt_i / dp_i
                        if freq_i < min_freq:
                            continue
                        pval_i = binom_upper_tail_pval(alt_i, dp_i, err_rate)
                        if pval_i <= p_cutoff and (not passed or pval_i < best_pval):
                            passed = True
                            best_dp = dp_i
                            best_alt = alt_i
                            best_freq = freq_i
                            best_pval = pval_i
                            best_source = bam_labels[bam_idx]

                if passed:
                    snp_records.append(
                        (
                            rec_chrom,
                            pos_int,
                            ref,
                            alt_base,
                            best_dp,
                            best_alt,
                            best_freq,
                            best_pval,
                            best_source,
                        )
                    )

        ret = proc.wait()
        if ret != 0:
            err_fh.seek(0)
            stderr_txt = err_fh.read()
            raise RuntimeError(
                f"samtools mpileup failed on region {region_str} (exit {ret}): {stderr_txt.strip()}"
            )

    return snp_records


def write_empty_outputs(output_txt, output_tsv):
    os.makedirs(os.path.dirname(os.path.abspath(output_txt)), exist_ok=True)
    with open(output_txt, "w", newline="\n", encoding="utf-8"):
        pass

    if output_tsv:
        os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True)
        with open(output_tsv, "w", newline="", encoding="utf-8") as tsv_fh:
            wr = csv.writer(tsv_fh, delimiter="\t", lineterminator="\n")
            wr.writerow(
                ["chrom", "pos", "ref", "alt", "depth", "alt_count", "freq", "pval", "source"]
            )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Call mutational anomalies / SNPs from BAM pileups using exact binomial & frequency filters"
    )
    parser.add_argument(
        "-b",
        "--bam",
        nargs="*",
        default=[],
        help="Coordinate-sorted, indexed BAM file(s) to inspect",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=False,
        default="",
        help="Reference genome FASTA file",
    )
    parser.add_argument(
        "--fai",
        required=False,
        default="",
        help="Reference genome FASTA .fai index (defaults to <fasta>.fai)",
    )
    parser.add_argument(
        "-o",
        "--output-txt",
        required=True,
        help="Output snp.txt path (REF:ALT:CHR:POS format)",
    )
    parser.add_argument(
        "--output-tsv",
        default="",
        help="Optional diagnostic TSV path with depth, frequency, and p-value columns",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of parallel samtools mpileup worker processes (default: 1)",
    )
    parser.add_argument(
        "--min-qual",
        type=int,
        default=62,
        help="fastq2EZbakR minqual threshold; converted to Phred >= max(0, min_qual - 32) (default: 62 -> Phred 30)",
    )
    parser.add_argument(
        "--min-bq",
        type=int,
        default=None,
        help="Optional direct Phred base quality threshold overriding --min-qual",
    )
    parser.add_argument(
        "--min-depth",
        type=int,
        default=5,
        help="Minimum non-intron base coverage required to evaluate a site (default: 5)",
    )
    parser.add_argument(
        "--min-alt",
        type=int,
        default=3,
        help="Minimum alternate base observations required (default: 3)",
    )
    parser.add_argument(
        "--min-freq",
        type=float,
        default=0.25,
        help="Minimum alternate allele frequency alt_count / depth (default: 0.25)",
    )
    parser.add_argument(
        "--err-rate",
        type=float,
        default=0.002,
        help="Null background error/mutation rate p0 for upper-tail binomial test (default: 0.002)",
    )
    parser.add_argument(
        "--p-cutoff",
        type=float,
        default=1e-4,
        help="Maximum upper-tail binomial p-value to flag an anomaly (default: 1e-4)",
    )
    parser.add_argument(
        "--mpileup-extra",
        nargs="?",
        const="",
        default="",
        type=str,
        help="Optional extra flags passed to samtools mpileup (use --mpileup-extra='...' if starting with '-')",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    bam_paths = [b for b in args.bam if b and str(b).strip()]

    print(f"Start call_snps: {datetime.datetime.now()}", flush=True)

    if not bam_paths:
        write_empty_outputs(args.output_txt, args.output_tsv)
        print("No input BAMs provided for anomaly/SNP calling; wrote empty snp.txt.", flush=True)
        return

    if not args.fasta:
        raise ValueError("--fasta is required when one or more --bam files are provided.")

    fai_path = args.fai if args.fai else f"{args.fasta}.fai"
    if not os.path.exists(fai_path):
        raise FileNotFoundError(f"FASTA index (.fai) not found: {fai_path}")

    min_phred = args.min_bq if args.min_bq is not None else max(0, args.min_qual - 32)
    min_depth = max(1, args.min_depth)
    min_alt = max(1, args.min_alt)
    min_freq = float(args.min_freq)
    err_rate = float(args.err_rate)
    p_cutoff = float(args.p_cutoff)
    threads = max(1, args.threads)
    extra_args = shlex.split(args.mpileup_extra) if args.mpileup_extra.strip() else []

    bam_labels = [
        os.path.splitext(os.path.basename(p))[0] for p in bam_paths
    ]

    all_contigs = read_fai_contigs(fai_path)
    active_contigs = get_active_contigs(all_contigs, bam_paths)
    regions = build_regions(active_contigs, threads)

    print(
        f"Inspecting {len(bam_paths)} BAM(s) across {len(active_contigs)}/{len(all_contigs)} active contigs "
        f"({len(regions)} region chunks) using {threads} worker(s) "
        f"[Phred>={min_phred}, depth>={min_depth}, alt>={min_alt}, freq>={min_freq}, p0={err_rate}, p<={p_cutoff}].",
        flush=True,
    )

    if not regions:
        write_empty_outputs(args.output_txt, args.output_tsv)
        print("Zero mapped reads found across input BAMs; wrote empty snp.txt.", flush=True)
        return

    tasks = [
        (
            reg,
            bam_paths,
            bam_labels,
            args.fasta,
            min_phred,
            min_depth,
            min_alt,
            min_freq,
            err_rate,
            p_cutoff,
            extra_args,
        )
        for reg in regions
    ]

    all_records = []
    n_workers = min(threads, len(tasks))
    if n_workers <= 1:
        for t in tasks:
            all_records.extend(_scan_region_worker(t))
    else:
        with mp.Pool(processes=n_workers) as pool:
            for chunk_records in pool.imap_unordered(_scan_region_worker, tasks, chunksize=1):
                all_records.extend(chunk_records)

    # Sort deterministically by (chrom, pos, ref, alt) and deduplicate snp.txt entries
    all_records.sort(key=lambda r: (r[0], r[1], r[2], r[3]))

    os.makedirs(os.path.dirname(os.path.abspath(args.output_txt)), exist_ok=True)
    snp_strings = sorted({f"{r[2]}:{r[3]}:{r[0]}:{r[1]}" for r in all_records})
    with open(args.output_txt, "w", newline="\n", encoding="utf-8") as txt_fh:
        for entry in snp_strings:
            txt_fh.write(f"{entry}\n")

    if args.output_tsv:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_tsv)), exist_ok=True)
        with open(args.output_tsv, "w", newline="", encoding="utf-8") as tsv_fh:
            wr = csv.writer(tsv_fh, delimiter="\t", lineterminator="\n")
            wr.writerow(
                ["chrom", "pos", "ref", "alt", "depth", "alt_count", "freq", "pval", "source"]
            )
            for rec in all_records:
                wr.writerow(
                    [
                        rec[0],
                        rec[1],
                        rec[2],
                        rec[3],
                        rec[4],
                        rec[5],
                        f"{rec[6]:.6f}",
                        f"{rec[7]:.4e}",
                        rec[8],
                    ]
                )

    print(
        f"Finished call_snps: {len(snp_strings)} anomalous alleles written to {args.output_txt} "
        f"({datetime.datetime.now()}).",
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR in call_snps.py: {exc}", file=sys.stderr, flush=True)
        raise
