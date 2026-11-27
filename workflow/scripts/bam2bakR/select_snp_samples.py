#!/usr/bin/env python3
"""
    File name: select_snp_samples.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Rank candidate filtered BAMs by mapped read count and select up to
                 max_samples top-coverage samples for no-control anomaly calling.
"""

import argparse
import concurrent.futures
import csv
import gzip
import os
import shutil
import struct
import subprocess
import sys


def _infer_sample_name(bam_path):
    """Extract sample ID from `results/sf_reads/{sample}.s.bam` or `{sample}.bam`."""
    base = os.path.basename(bam_path)
    if base.endswith(".s.bam"):
        return base[:-6]
    if base.endswith(".bam"):
        return base[:-4]
    return os.path.splitext(base)[0]


def _read_exact(handle, nbytes):
    """Read up to `nbytes` uncompressed bytes from a binary stream."""
    chunks = []
    remaining = nbytes
    while remaining > 0:
        chunk = handle.read(remaining)
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def _skip_exact(handle, nbytes):
    """Advance `nbytes` uncompressed bytes or raise ValueError on premature EOF."""
    remaining = nbytes
    while remaining > 0:
        chunk = handle.read(min(65536, remaining))
        if not chunk:
            raise ValueError("Unexpected EOF while scanning BAM record stream.")
        remaining -= len(chunk)


def _count_bam_records_python(bam_path):
    """Pure-Python BGZF/BAM record counter used when `samtools` is unavailable."""
    with gzip.open(bam_path, "rb") as handle:
        magic = _read_exact(handle, 4)
        if magic != b"BAM\x01":
            raise ValueError(f"Not a valid BAM file: {bam_path}")

        l_text_raw = _read_exact(handle, 4)
        if len(l_text_raw) != 4:
            raise ValueError(f"Truncated BAM header in {bam_path}")
        l_text = struct.unpack("<I", l_text_raw)[0]
        _skip_exact(handle, l_text)

        n_ref_raw = _read_exact(handle, 4)
        if len(n_ref_raw) != 4:
            raise ValueError(f"Truncated BAM reference table in {bam_path}")
        n_ref = struct.unpack("<I", n_ref_raw)[0]

        for _ in range(n_ref):
            l_name_raw = _read_exact(handle, 4)
            if len(l_name_raw) != 4:
                raise ValueError(f"Truncated BAM reference entry in {bam_path}")
            l_name = struct.unpack("<I", l_name_raw)[0]
            _skip_exact(handle, l_name + 4)

        count = 0
        while True:
            rec_hdr = _read_exact(handle, 4)
            if not rec_hdr:
                break
            if len(rec_hdr) != 4:
                raise ValueError(f"Truncated BAM alignment record in {bam_path}")
            block_size = struct.unpack("<I", rec_hdr)[0]
            _skip_exact(handle, block_size)
            count += 1

    return count


def count_bam_records(bam_path, extra_threads=0, use_samtools=True):
    """Count alignment records in `bam_path` via `samtools view -c` or Python."""
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    if use_samtools and shutil.which("samtools"):
        cmd = ["samtools", "view", "-c"]
        if extra_threads > 0:
            cmd.extend(["-@", str(extra_threads)])
        cmd.append(bam_path)
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"samtools view -c failed on {bam_path} "
                f"(exit {proc.returncode}): {proc.stderr.strip()}"
            )
        return int(proc.stdout.strip())

    return _count_bam_records_python(bam_path)


def rank_and_select_samples(
    bam_paths,
    max_samples,
    samples=None,
    threads=1,
    use_samtools=True,
):
    """
    Count records in each BAM and return `(selected_samples, ranked_rows)`, where
    `ranked_rows` is sorted by `(-mapped_records, sample)`.
    """
    if not bam_paths:
        raise ValueError("At least one BAM path must be provided.")
    if isinstance(max_samples, bool) or (max_samples != -1 and max_samples < 1):
        raise ValueError(
            f"Invalid max_samples '{max_samples}'. Expected -1 or an integer >= 1."
        )

    if samples is None:
        sample_names = [_infer_sample_name(p) for p in bam_paths]
    else:
        if len(samples) != len(bam_paths):
            raise ValueError(
                f"Number of sample names ({len(samples)}) does not match "
                f"number of BAM paths ({len(bam_paths)})."
            )
        sample_names = [str(s).strip() for s in samples]

    if len(set(sample_names)) != len(sample_names):
        raise ValueError(f"Duplicate sample names provided: {sample_names}")

    total_threads = max(1, int(threads))
    n_workers = min(total_threads, len(bam_paths))
    extra_threads = max(0, (total_threads // n_workers) - 1)

    def _worker(item):
        sample_id, path = item
        depth = count_bam_records(
            path, extra_threads=extra_threads, use_samtools=use_samtools
        )
        return sample_id, depth

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as pool:
        counts = list(pool.map(_worker, zip(sample_names, bam_paths)))

    ranked = sorted(counts, key=lambda pair: (-pair[1], pair[0]))
    cap = len(ranked) if max_samples == -1 else min(max_samples, len(ranked))
    selected = [sample_id for sample_id, _ in ranked[:cap]]
    selected_set = set(selected)
    ranked_rows = [
        (sample_id, depth, sample_id in selected_set) for sample_id, depth in ranked
    ]
    return selected, ranked_rows


def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Select top-coverage samples from filtered BAMs for no-control "
            "mutational-anomaly calling."
        )
    )
    parser.add_argument(
        "-b",
        "--bams",
        nargs="+",
        required=True,
        help="Filtered BAM file(s) (`results/sf_reads/{sample}.s.bam`).",
    )
    parser.add_argument(
        "-s",
        "--samples",
        nargs="+",
        default=None,
        help="Optional explicit sample names matching --bams order.",
    )
    parser.add_argument(
        "-m",
        "--max-samples",
        type=int,
        required=True,
        help="Maximum number of samples to select (-1 = all, or integer >= 1).",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Total thread budget across parallel BAM record counts (default: 1).",
    )
    parser.add_argument(
        "-o",
        "--output-selected",
        required=True,
        help="Output newline-delimited file of selected sample names.",
    )
    parser.add_argument(
        "--output-depths",
        default="",
        help="Optional output TSV recording sample, mapped_records, and selected.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    selected, ranked_rows = rank_and_select_samples(
        args.bams,
        args.max_samples,
        samples=args.samples,
        threads=args.threads,
    )

    sel_dir = os.path.dirname(os.path.abspath(args.output_selected))
    if sel_dir:
        os.makedirs(sel_dir, exist_ok=True)
    with open(args.output_selected, "w", encoding="utf-8", newline="\n") as fh:
        for sample_id in selected:
            fh.write(f"{sample_id}\n")

    if args.output_depths:
        dep_dir = os.path.dirname(os.path.abspath(args.output_depths))
        if dep_dir:
            os.makedirs(dep_dir, exist_ok=True)
        with open(args.output_depths, "w", encoding="utf-8", newline="") as tsv_fh:
            writer = csv.writer(tsv_fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", "mapped_records", "selected"])
            for sample_id, depth, is_sel in ranked_rows:
                writer.writerow([sample_id, depth, is_sel])

    chosen_summary = ", ".join(
        f"{sample_id} ({depth})"
        for sample_id, depth, is_sel in ranked_rows
        if is_sel
    )
    print(
        f"Selected {len(selected)}/{len(ranked_rows)} sample(s) by filtered read "
        f"depth: {chosen_summary}",
        flush=True,
    )


if __name__ == "__main__":
    main()
