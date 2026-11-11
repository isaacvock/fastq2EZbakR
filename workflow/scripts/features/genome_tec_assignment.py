#!/usr/bin/env python3
"""
    File name: genome_tec_assignment.py
    Author: Isaac Vock
    Email: isaac.vock@gmail.com
    License: GPLv3
    Description: Aligner-agnostic, strand- and splice-aware assignment of genome-aligned
                 reads/read-pairs to transcript equivalence classes (TECs).
"""

import argparse
from bisect import bisect_right
import csv
import gzip
import multiprocessing as mp
import os
import shutil
import struct
import tempfile
import zlib

import pysam


BIN_SHIFT = 14  # 16,384 bp spatial index bins
MIN_PARALLEL_BAM_BYTES = 2 * 1024 * 1024

_GLOBAL_TX_INDEX = None


class TranscriptModel:
    __slots__ = (
        "transcript_id",
        "exon_starts",
        "exon_ends",
        "introns",
        "eff_tx_start",
        "eff_tx_end",
    )

    def __init__(self, transcript_id, exon_starts, exon_ends, introns, end_overhang):
        self.transcript_id = transcript_id
        self.exon_starts = exon_starts
        self.exon_ends = exon_ends
        self.introns = introns
        self.eff_tx_start = exon_starts[0] - end_overhang
        self.eff_tx_end = exon_ends[-1] + end_overhang


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


def load_transcript_index(gtf_path, end_overhang=2):
    """
    Parse `exon` records from `gtf_path` into binned per-`(chrom, strand)` transcript models.
    Coordinates are stored as 0-based half-open `[start, end)`.
    """
    opener = (
        (lambda p: gzip.open(p, "rt", encoding="utf-8"))
        if str(gtf_path).endswith(".gz")
        else (lambda p: open(p, "r", encoding="utf-8"))
    )

    raw_exons = {}
    with opener(gtf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            strand = fields[6]
            if strand not in ("+", "-"):
                continue
            tx_id = _extract_transcript_id(fields[8])
            if not tx_id:
                continue
            start_0 = int(fields[3]) - 1
            end_0 = int(fields[4])
            if end_0 <= start_0:
                continue

            key = (fields[0], strand, tx_id)
            ex_list = raw_exons.get(key)
            if ex_list is None:
                raw_exons[key] = [(start_0, end_0)]
            else:
                ex_list.append((start_0, end_0))

    bins_by_cs = {}
    for (chrom, strand, tx_id), ex_list in raw_exons.items():
        ex_list.sort()
        merged = [ex_list[0]]
        for s, e in ex_list[1:]:
            prev_s, prev_e = merged[-1]
            if s <= prev_e:
                if e > prev_e:
                    merged[-1] = (prev_s, e)
            else:
                merged.append((s, e))

        exon_starts = tuple(s for s, _ in merged)
        exon_ends = tuple(e for _, e in merged)
        introns = frozenset(
            (exon_ends[i], exon_starts[i + 1]) for i in range(len(merged) - 1)
        )
        model = TranscriptModel(tx_id, exon_starts, exon_ends, introns, end_overhang)

        cs_bins = bins_by_cs.setdefault((chrom, strand), {})
        b_start = max(0, model.eff_tx_start) >> BIN_SHIFT
        b_end = max(0, model.eff_tx_end) >> BIN_SHIFT
        for b in range(b_start, b_end + 1):
            bucket = cs_bins.get(b)
            if bucket is None:
                cs_bins[b] = [model]
            else:
                bucket.append(model)

    return bins_by_cs


def _parse_bgzf_block_at(data, offset):
    n = len(data)
    if offset + 28 > n:
        return None
    if (
        data[offset] != 0x1F
        or data[offset + 1] != 0x8B
        or data[offset + 2] != 0x08
    ):
        return None
    flg = data[offset + 3]
    if (flg & 0x04) == 0 or (flg & 0xE0) != 0:
        return None

    xlen = struct.unpack_from("<H", data, offset + 10)[0]
    extra_end = offset + 12 + xlen
    if extra_end + 8 > n:
        return None

    bsize = None
    sub_pos = offset + 12
    while sub_pos + 4 <= extra_end:
        si1, si2, slen = struct.unpack_from("<BBH", data, sub_pos)
        if si1 == 66 and si2 == 67 and slen == 2 and sub_pos + 6 <= extra_end:
            bsize = struct.unpack_from("<H", data, sub_pos + 4)[0]
        sub_pos += 4 + slen

    if bsize is None:
        return None

    block_len = bsize + 1
    block_end = offset + block_len
    if block_end > n or block_len < 28:
        return None

    payload_start = extra_end
    if flg & 0x08:
        zero_idx = data.find(b"\x00", payload_start, block_end - 8)
        if zero_idx < 0:
            return None
        payload_start = zero_idx + 1
    if flg & 0x10:
        zero_idx = data.find(b"\x00", payload_start, block_end - 8)
        if zero_idx < 0:
            return None
        payload_start = zero_idx + 1
    if flg & 0x02:
        payload_start += 2

    if payload_start > block_end - 8:
        return None

    crc_expected, isize_expected = struct.unpack_from("<II", data, block_end - 8)
    if isize_expected > 65536:
        return None

    try:
        raw = zlib.decompress(data[payload_start : block_end - 8], -15)
    except zlib.error:
        return None

    if len(raw) != isize_expected:
        return None
    if (zlib.crc32(raw) & 0xFFFFFFFF) != crc_expected:
        return None

    return block_len, raw


def _is_valid_bam_record(ubuf, pos, ref_lengths):
    n = len(ubuf)
    if pos + 36 > n:
        return "TRUNCATED"

    (
        block_size,
        ref_id,
        aln_pos,
        l_read_name,
        _mapq,
        _bin,
        n_cigar_op,
        flag,
        l_seq,
        next_ref_id,
        next_pos,
        _tlen,
    ) = struct.unpack_from("<iiiBBHHHiiii", ubuf, pos)

    n_refs = len(ref_lengths)
    if l_read_name < 2 or l_seq < 0 or (flag & 0xF000) != 0:
        return None
    if ref_id < -1 or ref_id >= n_refs or next_ref_id < -1 or next_ref_id >= n_refs:
        return None
    if aln_pos < -1 or (ref_id >= 0 and aln_pos >= ref_lengths[ref_id]):
        return None
    if next_pos < -1 or (next_ref_id >= 0 and next_pos >= ref_lengths[next_ref_id]):
        return None

    min_body = 32 + l_read_name + 4 * n_cigar_op + ((l_seq + 1) // 2) + l_seq
    if block_size < min_body or block_size > 64 * 1024 * 1024:
        return None

    rec_end = pos + 4 + block_size
    cigar_end = pos + 36 + l_read_name + 4 * n_cigar_op
    if cigar_end > n:
        return "TRUNCATED"

    rname_bytes = ubuf[pos + 36 : pos + 36 + l_read_name]
    if rname_bytes[-1] != 0:
        return None
    for b in rname_bytes[:-1]:
        if b < 33 or b > 126:
            return None

    if n_cigar_op > 0:
        q_sum = 0
        cigar_offset = pos + 36 + l_read_name
        for idx in range(n_cigar_op):
            cigar_int = struct.unpack_from("<I", ubuf, cigar_offset + 4 * idx)[0]
            op = cigar_int & 0xF
            oplen = cigar_int >> 4
            if op > 8 or oplen <= 0:
                return None
            if op in (0, 1, 4, 7, 8):
                q_sum += oplen
        if l_seq > 0 and q_sum != l_seq:
            return None

    if rec_end > n:
        return "TRUNCATED"
    return rec_end


def _find_voffset_near(fh, file_size, min_block_offset, target_byte, ref_lengths):
    search_pos = max(min_block_offset, target_byte)
    if search_pos >= file_size - 28:
        return None

    fh.seek(search_pos)
    window = fh.read(min(512 * 1024, file_size - search_pos))

    b0_rel = None
    scan = 0
    limit = len(window) - 28
    while scan <= limit:
        idx = window.find(b"\x1f\x8b\x08", scan)
        if idx < 0 or idx > limit:
            break
        parsed0 = _parse_bgzf_block_at(window, idx)
        if parsed0 is not None:
            blen0, _ = parsed0
            next1 = idx + blen0
            if search_pos + next1 == file_size:
                b0_rel = idx
                break
            parsed1 = _parse_bgzf_block_at(window, next1)
            if parsed1 is not None:
                blen1, _ = parsed1
                next2 = next1 + blen1
                if (
                    search_pos + next2 == file_size
                    or _parse_bgzf_block_at(window, next2) is not None
                ):
                    b0_rel = idx
                    break
        scan = idx + 1

    if b0_rel is None:
        return None

    cur_file_offset = search_pos + b0_rel
    while cur_file_offset < file_size - 28:
        fh.seek(cur_file_offset)
        buf = fh.read(min(512 * 1024, file_size - cur_file_offset))
        blocks = []
        rel = 0
        reached_eof = False
        while rel < len(buf) and len(blocks) < 6:
            parsed = _parse_bgzf_block_at(buf, rel)
            if parsed is None:
                break
            blen, raw = parsed
            if len(raw) == 0 and cur_file_offset + rel + blen == file_size:
                reached_eof = True
                break
            blocks.append((cur_file_offset + rel, blen, raw))
            rel += blen
            if cur_file_offset + rel == file_size:
                reached_eof = True
                break

        if not blocks:
            return None

        first_block_offset, first_blen, first_raw = blocks[0]
        if len(first_raw) == 0:
            cur_file_offset = first_block_offset + first_blen
            continue

        ubuf = b"".join(b[2] for b in blocks)
        first_u_len = len(first_raw)

        for uoff in range(first_u_len):
            cur_p = uoff
            valid_chain = 0
            ok = True
            while valid_chain < 4:
                nxt = _is_valid_bam_record(ubuf, cur_p, ref_lengths)
                if nxt is None:
                    ok = False
                    break
                if nxt == "TRUNCATED":
                    ok = valid_chain >= 2
                    break
                valid_chain += 1
                cur_p = nxt
                if cur_p == len(ubuf) and reached_eof:
                    break
            if ok and valid_chain >= 1:
                return (first_block_offset << 16) | uoff

        cur_file_offset = first_block_offset + first_blen

    return None


def compute_bam_chunks(bam_path, threads):
    with pysam.AlignmentFile(bam_path, "rb") as sf:
        first_voff = sf.tell()
        try:
            next(sf)
        except StopIteration:
            return [], True
        ref_lengths = list(sf.lengths)

    file_size = os.path.getsize(bam_path)
    if threads <= 1 or file_size < MIN_PARALLEL_BAM_BYTES:
        return [(first_voff, None)], False

    min_block_offset = (first_voff >> 16) + 1
    split_voffs = [first_voff]

    with open(bam_path, "rb") as fh, pysam.AlignmentFile(bam_path, "rb") as verify_sf:
        for k in range(1, threads):
            target = (k * file_size) // threads
            voff = _find_voffset_near(fh, file_size, min_block_offset, target, ref_lengths)
            if voff is None or voff <= split_voffs[-1]:
                continue
            try:
                verify_sf.seek(voff)
                rec = next(verify_sf)
                if not rec.query_name:
                    continue
            except (OSError, StopIteration, ValueError):
                continue
            split_voffs.append(voff)

    chunks = []
    for idx, start_voff in enumerate(split_voffs):
        end_voff = split_voffs[idx + 1] if idx + 1 < len(split_voffs) else None
        chunks.append((start_voff, end_voff))
    return chunks, False


def _infer_rna_strand(read, norm_strand):
    aligned_strand = "-" if read.is_reverse else "+"
    if read.is_paired:
        if read.is_read1:
            same = norm_strand == "forward"
        elif read.is_read2:
            same = norm_strand == "reverse"
        else:
            return None
    else:
        same = norm_strand == "forward"
    return aligned_strand if same else ("-" if aligned_strand == "+" else "+")


def _extract_blocks_and_introns(read, out_blocks, out_introns):
    """
    Extract 0-based half-open exon alignment blocks `(start, end)` and CIGAR `N`
    intron intervals `(start, end)` from `read`.
    Merges intra-exonic insertions/deletions (`I`, `D`) into contiguous exon blocks.
    """
    cigartuples = read.cigartuples
    if not cigartuples:
        return

    r_pos = read.reference_start
    block_start = None

    for op, length in cigartuples:
        if op in (0, 2, 7, 8):  # M, D, =, X consume exonic reference
            if block_start is None:
                block_start = r_pos
            r_pos += length
        elif op == 3:  # N splice skip closes exon block and records intron
            if block_start is not None and r_pos > block_start:
                out_blocks.append((block_start, r_pos))
                block_start = None
            out_introns.add((r_pos, r_pos + length))
            r_pos += length
        # op in (1=I, 4=S, 5=H, 6=P): no reference consumption

    if block_start is not None and r_pos > block_start:
        out_blocks.append((block_start, r_pos))


def _blocks_fit_transcript(blocks, model, end_overhang):
    exon_starts = model.exon_starts
    exon_ends = model.exon_ends
    last_idx = len(exon_starts) - 1

    for b_start, b_end in blocks:
        idx = bisect_right(exon_starts, b_start) - 1
        if idx < 0:
            idx = 0
        eff_s = exon_starts[idx] - (end_overhang if idx == 0 else 0)
        eff_e = exon_ends[idx] + (end_overhang if idx == last_idx else 0)
        if b_start < eff_s or b_end > eff_e:
            return False
    return True


def _assign_group_tecs(valid_reads, tx_index, norm_strand, end_overhang):
    chrom = valid_reads[0].reference_name
    for r in valid_reads[1:]:
        if r.reference_name != chrom:
            return None

    if norm_strand == "unstranded":
        candidate_strands = ("+", "-")
    else:
        inferred = _infer_rna_strand(valid_reads[0], norm_strand)
        if inferred is None:
            return None
        for r in valid_reads[1:]:
            if _infer_rna_strand(r, norm_strand) != inferred:
                return None
        candidate_strands = (inferred,)

    blocks = []
    introns = set()
    for r in valid_reads:
        _extract_blocks_and_introns(r, blocks, introns)

    if not blocks:
        return None

    min_pos = min(b[0] for b in blocks)
    max_pos = max(b[1] for b in blocks)
    bin_min = max(0, min_pos) >> BIN_SHIFT
    bin_max = max(0, max_pos) >> BIN_SHIFT

    matched_tx = set()
    for strand in candidate_strands:
        cs_bins = tx_index.get((chrom, strand))
        if not cs_bins:
            continue
        b_list_min = cs_bins.get(bin_min)
        if not b_list_min:
            continue
        b_list_max = cs_bins.get(bin_max)
        if not b_list_max:
            continue

        candidates = b_list_min if len(b_list_min) <= len(b_list_max) else b_list_max
        for model in candidates:
            if min_pos < model.eff_tx_start or max_pos > model.eff_tx_end:
                continue
            if introns and not introns.issubset(model.introns):
                continue
            if _blocks_fit_transcript(blocks, model, end_overhang):
                matched_tx.add(model.transcript_id)

    return matched_tx if matched_tx else None


def _init_pool_worker(gtf_path, end_overhang):
    global _GLOBAL_TX_INDEX
    if _GLOBAL_TX_INDEX is None:
        _GLOBAL_TX_INDEX = load_transcript_index(gtf_path, end_overhang=end_overhang)


def _process_chunk_worker(task):
    (
        worker_idx,
        bam_path,
        start_voff,
        end_voff,
        part_path,
        norm_strand,
        end_overhang,
    ) = task

    tx_index = _GLOBAL_TX_INDEX
    assert tx_index is not None

    with pysam.AlignmentFile(bam_path, "rb") as sf, open(
        part_path, "w", newline="\n", encoding="utf-8"
    ) as out_fh:
        wr = csv.writer(out_fh, lineterminator="\n")
        sf.seek(start_voff)

        first_rec = None
        first_rec_voff = None
        if worker_idx > 0:
            try:
                r0 = next(sf)
            except StopIteration:
                return
            skip_qname = r0.query_name
            while True:
                rec_voff = sf.tell()
                try:
                    r = next(sf)
                except StopIteration:
                    return
                if r.query_name != skip_qname:
                    if end_voff is not None and rec_voff > end_voff:
                        return
                    first_rec = r
                    first_rec_voff = rec_voff
                    break
        else:
            rec_voff = sf.tell()
            try:
                first_rec = next(sf)
                first_rec_voff = rec_voff
            except StopIteration:
                return

        cur_qname = first_rec.query_name
        cur_group = [first_rec]
        passed_end = end_voff is not None and first_rec_voff >= end_voff

        while True:
            rec_voff = sf.tell()
            try:
                rec = next(sf)
            except StopIteration:
                break

            qname = rec.query_name
            if qname == cur_qname:
                cur_group.append(rec)
                if end_voff is not None and rec_voff >= end_voff:
                    passed_end = True
            else:
                valid = [
                    r
                    for r in cur_group
                    if not r.is_unmapped
                    and not r.is_secondary
                    and not r.is_supplementary
                    and r.reference_name is not None
                ]
                if valid:
                    tecs = _assign_group_tecs(valid, tx_index, norm_strand, end_overhang)
                    if tecs:
                        wr.writerow([cur_qname, "+".join(sorted(tecs))])

                if passed_end or (end_voff is not None and rec_voff > end_voff):
                    cur_group = None
                    break

                cur_qname = qname
                cur_group = [rec]
                if end_voff is not None and rec_voff >= end_voff:
                    passed_end = True

        if cur_group:
            valid = [
                r
                for r in cur_group
                if not r.is_unmapped
                and not r.is_secondary
                and not r.is_supplementary
                and r.reference_name is not None
            ]
            if valid:
                tecs = _assign_group_tecs(valid, tx_index, norm_strand, end_overhang)
                if tecs:
                    wr.writerow([cur_qname, "+".join(sorted(tecs))])


def run_genome_tec_assignment(
    bam_path,
    gtf_path,
    output_csv,
    strandedness="reverse",
    threads=1,
    end_overhang=2,
):
    global _GLOBAL_TX_INDEX

    norm_strand = normalize_strandedness(strandedness)
    threads = max(1, int(threads))
    end_overhang = max(0, int(end_overhang))

    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    chunks, is_empty = compute_bam_chunks(bam_path, threads)

    if is_empty:
        with open(output_csv, "w", newline="\n", encoding="utf-8") as out_fh:
            csv.writer(out_fh, lineterminator="\n").writerow(
                ["qname", "bamfile_transcripts"]
            )
        return

    _GLOBAL_TX_INDEX = load_transcript_index(gtf_path, end_overhang=end_overhang)

    tmp_dir = tempfile.mkdtemp(
        prefix="tec_parts_", dir=os.path.dirname(os.path.abspath(output_csv))
    )
    try:
        part_paths = [
            os.path.join(tmp_dir, f"part_{idx:04d}.csv") for idx in range(len(chunks))
        ]
        tasks = [
            (
                idx,
                bam_path,
                start_voff,
                end_voff,
                part_paths[idx],
                norm_strand,
                end_overhang,
            )
            for idx, (start_voff, end_voff) in enumerate(chunks)
        ]

        if len(tasks) <= 1:
            _process_chunk_worker(tasks[0])
        else:
            start_method = mp.get_start_method(allow_none=True) or mp.get_start_method()
            if start_method == "fork":
                with mp.Pool(processes=len(tasks)) as pool:
                    pool.map(_process_chunk_worker, tasks)
            else:
                with mp.Pool(
                    processes=len(tasks),
                    initializer=_init_pool_worker,
                    initargs=(gtf_path, end_overhang),
                ) as pool:
                    pool.map(_process_chunk_worker, tasks)

        with open(output_csv, "wb") as out_bin:
            out_bin.write(b"qname,bamfile_transcripts\n")
            for p_path in part_paths:
                if os.path.exists(p_path):
                    with open(p_path, "rb") as part_bin:
                        shutil.copyfileobj(part_bin, out_bin, length=1024 * 1024)
    finally:
        _GLOBAL_TX_INDEX = None
        shutil.rmtree(tmp_dir, ignore_errors=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assign queryname-sorted genome alignments to strand- & splice-compatible TECs"
    )
    parser.add_argument(
        "-b", "--bam", required=True, help="Queryname-sorted genome BAM file (.s.bam)"
    )
    parser.add_argument("-g", "--gtf", required=True, help="Annotation GTF file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV path")
    parser.add_argument(
        "-s",
        "--strandedness",
        default="reverse",
        help="Library strandedness: 'reverse', 'yes'/'forward', or 'no' (default: reverse)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Number of parallel workers"
    )
    parser.add_argument(
        "--end-overhang",
        type=int,
        default=2,
        help="Permitted overhang (bp) at annotated transcript 5'/3' ends (default: 2; internal splice sites always 0)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    cli = parse_args()
    run_genome_tec_assignment(
        cli.bam,
        cli.gtf,
        cli.output,
        strandedness=cli.strandedness,
        threads=cli.threads,
        end_overhang=cli.end_overhang,
    )
