"""Trim only cDNA and copy the corresponding original barcode FASTQ records."""

import gzip
import itertools
import shlex
import subprocess
import tempfile
from pathlib import Path


def records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as handle:
        number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            number += 1
            record = (header, handle.readline(), handle.readline(), handle.readline())
            seq, qual = record[1].rstrip(b"\r\n"), record[3].rstrip(b"\r\n")
            if (not header.startswith(b"@") or not header[1:].strip()
                    or not record[2].startswith(b"+") or not record[3]
                    or not seq or len(seq) != len(qual)):
                raise ValueError(f"Malformed FASTQ record {number} in {path}")
            yield record


def read_name(record):
    name = record[0][1:].split()[0]
    return name[:-2] if name.endswith((b"/1", b"/2")) else name


def pairs(paths):
    for number, pair in enumerate(itertools.zip_longest(*(records(p) for p in paths)), 1):
        if None in pair:
            raise ValueError(f"FASTQs have different record counts at pair {number}")
        if read_name(pair[0]) != read_name(pair[1]):
            raise ValueError(f"FASTQ names disagree at pair {number}: "
                             f"{pair[0][0].strip()!r} / {pair[1][0].strip()!r}")
        yield number, pair


def fastp_options(adapters, extra):
    options = shlex.split(adapters) + shlex.split(extra)
    # The wrapper owns identity, layout, I/O, and worker count. R1/SE options
    # apply to the cDNA, even when it came from the original R2 FASTQ.
    forbidden = {
        "--in1", "--in2", "--out1", "--out2", "--stdin", "--stdout",
        "--interleaved_in", "--unpaired1", "--unpaired2", "--failed_out",
        "--json", "--html", "--thread", "--split", "--split_by_lines",
        "--split_prefix_digits", "--reads_to_process", "--dedup",
        "--umi", "--merge", "--merged_out", "--include_unmerged",
        "--correction", "--fix_mgi_id", "--detect_adapter_for_pe",
        "--adapter_sequence_r2", "--trim_front2", "--trim_tail2", "--max_len2",
    }
    forbidden_short = set("iIoOhjwUsSmDcFTh2")
    for token in options:
        flag = token.split("=", 1)[0]
        if (flag in forbidden or flag.startswith("--umi_")
                or (flag.startswith("-") and not flag.startswith("--")
                    and len(flag) > 1 and flag[1] in forbidden_short)):
            raise ValueError(f"fastp option {flag} is incompatible with single-cell "
                             "preprocessing; use single-end cDNA trimming options only.")
    return options


def preprocess(paths, outputs, cdna_read, skip, adapters="", extra="",
               threads=1, html=None, json=None, log=None):
    if len(paths) != 2 or len(outputs) != 2 or cdna_read not in (1, 2):
        raise ValueError("Single-cell preprocessing requires two FASTQs and sc_cdna_read: 1 or 2")
    cdna = cdna_read - 1
    barcode = 1 - cdna
    if skip:
        with open(outputs[0], "wb") as out1, open(outputs[1], "wb") as out2:
            for _, pair in pairs(paths):
                out1.writelines(pair[0])
                out2.writelines(pair[1])
        return

    options = fastp_options(adapters, extra)
    marker = b":f2e_sc_"
    with tempfile.TemporaryDirectory(prefix="fastq2ezbakr-sc-") as scratch:
        staged = Path(scratch) / "cdna.fastq"
        trimmed = Path(scratch) / "trimmed.fastq"
        with staged.open("wb") as handle:
            for number, pair in pairs(paths):
                record = pair[cdna]
                name = record[0][1:].split(maxsplit=1)
                title = name[0] + marker + str(number).encode()
                if len(name) == 2:
                    title += b" " + name[1].rstrip(b"\r\n")
                handle.writelines((b"@" + title + b"\n", record[1], b"+\n", record[3]))
        subprocess.run([
            "fastp", "--in1", str(staged), "--out1", str(trimmed),
            "--thread", str(threads), "--html", str(html), "--json", str(json),
            *options,
        ], check=True, stdout=log, stderr=log)

        # fastp 0.23.4 preserves input order for file output, including with
        # multiple workers. Ordinals avoid matching a survivor to a different
        # occurrence of the same read name; no in-memory read-name set is needed.
        survivors = iter(records(trimmed))
        current = next(survivors, None)
        previous = 0

        def ordinal(record):
            return int(record[0][1:].split()[0].rsplit(marker, 1)[1])

        with open(outputs[cdna], "wb") as cdna_out, open(outputs[barcode], "wb") as barcode_out:
            for number, pair in pairs(paths):
                if current is None:
                    continue
                wanted = ordinal(current)
                if wanted <= previous or wanted < number:
                    raise ValueError("fastp changed or reordered single-cell read identifiers")
                if wanted == number:
                    original = pair[cdna]
                    expected = original[0][1:].split()[0] + marker + str(number).encode()
                    if current[0][1:].split()[0] != expected:
                        raise ValueError("fastp changed a single-cell read identifier")
                    cdna_out.writelines((original[0], current[1], original[2], current[3]))
                    barcode_out.writelines(pair[barcode])
                    previous = wanted
                    current = next(survivors, None)
            if current is not None:
                raise ValueError("fastp output contains a read absent from the input pair")


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log:
        preprocess(list(snakemake.input.fastqs),
                   [snakemake.output.r1, snakemake.output.r2],
                   snakemake.params.cdna_read, snakemake.params.skip,
                   snakemake.params.get("adapters", ""), snakemake.params.get("extra", ""),
                   snakemake.threads, snakemake.output.get("html"),
                   snakemake.output.get("json"), log)
