#!/usr/bin/env python3
"""
    File name: mut_call.py
    Author: Martin Machyna (parallel BGZF streaming, indel-flank masking & SNP denominator updates for fastq2EZbakR)
    Email: machyna@gmail.com
    Orcid: 0000-0002-3624-3472
    Wikidata: Q55137016
    License: GPLv3
    Description: Detect nucleotide recoding mutations directly from a queryname-sorted BAM file.
"""

import argparse
import csv
import datetime
import gzip
import heapq
import multiprocessing as mp
import os
import shutil
import struct
import tempfile
import zlib

import pysam


MUT_KEYS = (
    "TA", "CA", "GA", "NA",
    "AT", "CT", "GT", "NT",
    "AC", "TC", "GC", "NC",
    "AG", "TG", "CG", "NG",
    "AN", "TN", "CN", "GN", "NN",
)

DNA_COMP = {
    "A": "T", "C": "G", "T": "A", "G": "C", "N": "N",
    "a": "t", "c": "g", "t": "a", "g": "c", "n": "n",
}
VALID_REF_BASES = set(DNA_COMP.keys())
COMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CU_SPILL_LIMIT = 500_000
MIN_PARALLEL_BAM_BYTES = 2 * 1024 * 1024


def load_snps(snp_path):
    """Load called SNPs into a dict mapping chrom -> set of 0-based positions."""
    snp_by_chrom = {}
    if not snp_path or not os.path.exists(snp_path):
        return snp_by_chrom

    with open(snp_path, "r") as snp_file:
        for raw_line in snp_file:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split(":")
            if len(parts) < 4:
                continue
            chrom = ":".join(parts[2:-1])
            pos_0 = int(parts[-1]) - 1
            chrom_set = snp_by_chrom.get(chrom)
            if chrom_set is None:
                chrom_set = set()
                snp_by_chrom[chrom] = chrom_set
            chrom_set.add(pos_0)

    return snp_by_chrom


def _parse_bgzf_block_at(data, offset):
    """
    If `data[offset:]` starts with a valid BGZF block, return `(block_len, uncompressed_bytes)`.
    Otherwise return `None`.
    """
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
    if flg & 0x08:  # FNAME
        zero_idx = data.find(b"\x00", payload_start, block_end - 8)
        if zero_idx < 0:
            return None
        payload_start = zero_idx + 1
    if flg & 0x10:  # FCOMMENT
        zero_idx = data.find(b"\x00", payload_start, block_end - 8)
        if zero_idx < 0:
            return None
        payload_start = zero_idx + 1
    if flg & 0x02:  # FHCRC
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
    """
    Validate whether `ubuf[pos:]` begins a well-formed binary BAM alignment record.
    Returns `next_pos` (`pos + 4 + block_size`) if valid, `"TRUNCATED"` if valid header
    exceeds `len(ubuf)`, or `None` if invalid.
    """
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

    # Read name must be printable non-whitespace ASCII terminated by NUL
    rname_bytes = ubuf[pos + 36 : pos + 36 + l_read_name]
    if rname_bytes[-1] != 0:
        return None
    for b in rname_bytes[:-1]:
        if b < 33 or b > 126:
            return None

    # Validate CIGAR operations and query-consuming length
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
    """
    Find a verified BAM record virtual offset `(block_offset << 16) | uoffset` near `target_byte`.
    Returns `None` if no distinct valid interior record is found before EOF.
    """
    search_pos = max(min_block_offset, target_byte)
    if search_pos >= file_size - 28:
        return None

    fh.seek(search_pos)
    window = fh.read(min(512 * 1024, file_size - search_pos))

    # Locate the first valid BGZF block chain in `window`
    b0_rel = None
    scan = 0
    limit = len(window) - 28
    while scan <= limit:
        idx = window.find(b"\x1f\x8b\x08", scan)
        if idx < 0 or idx > limit:
            break
        parsed0 = _parse_bgzf_block_at(window, idx)
        if parsed0 is not None:
            blen0, raw0 = parsed0
            next1 = idx + blen0
            if search_pos + next1 == file_size:
                b0_rel = idx
                break
            parsed1 = _parse_bgzf_block_at(window, next1)
            if parsed1 is not None:
                blen1, _ = parsed1
                next2 = next1 + blen1
                if search_pos + next2 == file_size or _parse_bgzf_block_at(window, next2) is not None:
                    b0_rel = idx
                    break
        scan = idx + 1

    if b0_rel is None:
        return None

    # Decompress consecutive BGZF blocks starting at `cur_file_offset` until we find a record boundary
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

        # Entire first block was interior to an extra-long record; step to the next BGZF block
        cur_file_offset = first_block_offset + first_blen

    return None


def compute_bam_chunks(bam_path, threads):
    """
    Compute verified `[(start_voff, end_voff), ...]` chunks across `bam_path`.
    Returns `([], True)` if the BAM contains zero alignment records.
    """
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


def _flush_cu_spill(cU, spill_dir, worker_idx, spill_idx):
    """Sort `cU` in memory and spill to a compact temporary CSV file."""
    spill_path = os.path.join(spill_dir, f"cU_{worker_idx}_{spill_idx}.tsv")
    with open(spill_path, "w", newline="") as fh:
        for (rname, gloc, mt), (trials, n) in sorted(cU.items()):
            fh.write(f"{rname}\t{gloc}\t{mt}\t{trials}\t{n}\n")
    cU.clear()
    return spill_path


def _iter_cu_spill(spill_path):
    """Yield `((rname, gloc, mt), trials, n)` from a sorted spill file."""
    with open(spill_path, "r") as fh:
        for line in fh:
            rname, gloc_s, mt, trials_s, n_s = line.rstrip("\n").split("\t")
            yield (rname, int(gloc_s), mt), int(trials_s), int(n_s)


def _merge_cu_spills(spill_paths, cu_output_path):
    """Stream-merge sorted `cU` spill files into the final `_cU.csv` or `_cU.csv.gz`."""
    is_gz = cu_output_path.endswith(".gz")
    opener = (lambda p: gzip.open(p, "wt", newline="", compresslevel=3)) if is_gz else (lambda p: open(p, "w", newline=""))

    streams = [_iter_cu_spill(p) for p in spill_paths]
    merged = heapq.merge(*streams, key=lambda item: item[0])

    with opener(cu_output_path) as out_fh:
        wr = csv.writer(out_fh)
        wr.writerow(["rname", "gloc", "tp", "trials", "n"])

        cur_key = None
        cur_trials = 0
        cur_n = 0
        for key, trials, n in merged:
            if key == cur_key:
                cur_trials += trials
                cur_n += n
            else:
                if cur_key is not None:
                    wr.writerow([cur_key[0], cur_key[1], cur_key[2], cur_trials, cur_n])
                cur_key = key
                cur_trials = trials
                cur_n = n

        if cur_key is not None:
            wr.writerow([cur_key[0], cur_key[1], cur_key[2], cur_trials, cur_n])


def _process_read_group(
    records,
    strand_check,
    min_qual,
    min_dist,
    indel_flank,
    snp_by_chrom,
    mut_pos,
    muts_rds,
    tracks,
    mut_types,
    mut_type_set,
    base_set,
    base_to_muttypes,
    wr,
    cU,
    freq,
    track_files,
):
    """Process all alignments sharing one `query_name` and emit at most one row."""
    r0 = records[0]
    qname = r0.query_name
    rname = r0.reference_name
    if rname is None:
        return

    MDstore = {}
    indel_masked_rpos = set()
    has_sj = False
    fr_orient = ""

    for r in records:
        if r.is_unmapped or r.cigarstring is None or r.query_sequence is None:
            continue

        cigar_str = r.cigarstring
        if "N" in cigar_str:
            has_sj = True

        is_rev_orient = (
            (r.is_paired and (r.is_read1 == (r.is_reverse == strand_check)))
            or (not r.is_paired and (r.is_reverse == strand_check))
        )
        fr_orient = "R" if is_rev_orient else "F"

        bad_qpos = None
        if indel_flank >= 0 and ("I" in cigar_str or "D" in cigar_str):
            bad_qpos = set()
            q_cur = 0
            r_cur = r.reference_start
            for op, length in r.cigartuples:
                if op == 0 or op == 7 or op == 8:  # M, =, X
                    q_cur += length
                    r_cur += length
                elif op == 1:  # I
                    for qp in range(max(0, q_cur - indel_flank), q_cur + length + indel_flank):
                        bad_qpos.add(qp)
                    for rp in range(r_cur - indel_flank, r_cur + indel_flank):
                        indel_masked_rpos.add(rp)
                    q_cur += length
                elif op == 2:  # D
                    for qp in range(max(0, q_cur - indel_flank), q_cur + indel_flank):
                        bad_qpos.add(qp)
                    for rp in range(r_cur - indel_flank, r_cur + length + indel_flank):
                        indel_masked_rpos.add(rp)
                    r_cur += length
                elif op == 3:  # N
                    r_cur += length
                elif op == 4:  # S
                    q_cur += length

        try:
            aligned_pairs = r.get_aligned_pairs(matches_only=True, with_seq=True)
        except ValueError as exc:
            raise RuntimeError(
                f"Alignment '{qname}' is missing a valid MD tag required for mutation calling: {exc}"
            ) from exc

        seq = r.query_sequence
        if is_rev_orient:
            seq = seq.translate(COMP_TABLE)
        quals = r.query_qualities
        if quals is None:
            quals = [60] * len(seq)

        q_start = r.query_alignment_start
        q_len = r.query_alignment_length

        mate_md = {}
        for qpos, rpos, rbase in aligned_pairs:
            if bad_qpos is not None and qpos in bad_qpos:
                continue
            rel_q = qpos - q_start
            left_d = rel_q + 1
            right_d = q_len - rel_q
            end_dist = left_d if left_d <= right_d else right_d

            if is_rev_orient:
                ref_b = DNA_COMP.get(rbase, "N" if rbase.isupper() else "n")
            else:
                ref_b = rbase if rbase in VALID_REF_BASES else ("N" if rbase.isupper() else "n")

            mate_md[rpos] = (ref_b, seq[qpos], quals[qpos], end_dist)

        if not MDstore:
            MDstore = mate_md
        else:
            for pos, data in mate_md.items():
                prev = MDstore.get(pos)
                if prev is None:
                    MDstore[pos] = data
                else:
                    prev_ref, _, prev_q, _ = prev
                    cur_ref, _, cur_q, _ = data
                    prev_mut = prev_ref.islower()
                    cur_mut = cur_ref.islower()
                    if (
                        (prev_q < cur_q and prev_mut == cur_mut)
                        or (prev_q < cur_q and prev_mut and not cur_mut and prev_q + 33 <= min_qual)
                        or (cur_mut and not prev_mut and cur_q + 33 > min_qual)
                    ):
                        MDstore[pos] = data

    if indel_masked_rpos:
        for pos in indel_masked_rpos:
            MDstore.pop(pos, None)

    if not MDstore:
        return

    chrom_snps = snp_by_chrom.get(rname)
    nA = 0
    nC = 0
    nT = 0
    nG = 0
    muts = dict.fromkeys(MUT_KEYS, 0)
    gmutloc = [] if mut_pos else None
    tp = [] if mut_pos else None

    for pos, (ref_b, read_b, qual, end_dist) in MDstore.items():
        if chrom_snps is not None and pos in chrom_snps:
            continue

        qual_ok = (qual + 33) > min_qual
        ref_upper = ref_b.upper()
        is_mut = ref_b.islower()

        if qual_ok:
            if ref_upper == "T":
                nT += 1
            elif ref_upper == "A":
                nA += 1
            elif ref_upper == "C":
                nC += 1
            elif ref_upper == "G":
                nG += 1

        if muts_rds and ref_upper in base_set:
            fkey = f"{rname}:{pos}"
            entry = freq.get(fkey)
            if entry is None:
                entry = [1, 0]
                freq[fkey] = entry
            else:
                entry[0] += 1
            if is_mut and (ref_upper + read_b) in mut_type_set:
                entry[1] += 1

        if mut_pos and qual_ok and ref_upper in base_to_muttypes:
            for mt in base_to_muttypes[ref_upper]:
                ckey = (rname, pos + 1, mt)
                centry = cU.get(ckey)
                if centry is None:
                    cU[ckey] = [1, 0]
                else:
                    centry[0] += 1

        if is_mut and qual_ok and (end_dist + 33) > min_dist:
            mtype = ref_upper + read_b
            if mtype in muts:
                muts[mtype] += 1
            if mut_pos and mtype in mut_type_set:
                cU[(rname, pos + 1, mtype)][1] += 1
                gmutloc.append(str(pos))
                tp.append(mtype)

    sj_str = "TRUE" if has_sj else "FALSE"
    row = [qname, nA, nC, nT, nG, rname, fr_orient, sj_str]
    row.extend(muts[k] for k in MUT_KEYS)
    if mut_pos:
        row.append("|".join(gmutloc))
        row.append("|".join(tp))
    wr.writerow(row)

    if tracks and track_files is not None:
        for idx, mut in enumerate(mut_types):
            limit = 5 if muts[mut] > 5 else muts[mut]
            for c in range(0, limit + 1):
                track_files[c + idx * 6].write(qname + "\n")


# Global worker configuration populated before Pool fork / passed on spawn
_WORKER_CFG = {}


def _init_worker(cfg):
    global _WORKER_CFG
    _WORKER_CFG = cfg


def _run_chunk_worker(task):
    """
    Process BAM records in `[start_voff, end_voff]` with exact `query_name` group handoff:
    - If `worker_idx > 0`, skip all records sharing the `query_name` at `start_voff` (owned by previous worker).
    - Continue reading up through `rec_voff <= end_voff`, finishing the entire `query_name` group
      that spans or starts at `end_voff`.
    """
    worker_idx, start_voff, end_voff, out_chunk_path = task
    cfg = _WORKER_CFG

    bam_path = cfg["bam_path"]
    write_header = worker_idx == 0
    header = cfg["header"]
    is_gz = cfg["is_gz"]
    strand_check = cfg["strand_check"]
    min_qual = cfg["min_qual"]
    min_dist = cfg["min_dist"]
    indel_flank = cfg["indel_flank"]
    snp_by_chrom = cfg["snp_by_chrom"]
    mut_pos = cfg["mut_pos"]
    muts_rds = cfg["muts_rds"]
    tracks = cfg["tracks"]
    mut_types = cfg["mut_types"]
    mut_type_set = cfg["mut_type_set"]
    base_set = cfg["base_set"]
    base_to_muttypes = cfg["base_to_muttypes"]
    spill_dir = cfg["spill_dir"]
    track_prefix = cfg["track_prefix"]

    opener = (
        (lambda p: gzip.open(p, "wt", newline="", compresslevel=3))
        if is_gz
        else (lambda p: open(p, "w", newline=""))
    )

    cU = {} if mut_pos else None
    freq = {} if muts_rds else None
    cu_spills = []
    spill_idx = 0

    track_files = None
    if tracks:
        track_files = []
        for mt in mut_types:
            for i in range(6):
                track_files.append(open(f"{track_prefix}_{mt}_{i}_reads.txt", "w"))

    try:
        with pysam.AlignmentFile(bam_path, "rb") as sf, opener(out_chunk_path) as out_fh:
            wr = csv.writer(out_fh)
            if write_header:
                wr.writerow(header)

            sf.seek(start_voff)

            # Skip initial `query_name` group when `worker_idx > 0`
            first_rec = None
            first_rec_voff = None
            if worker_idx > 0:
                rec_voff = sf.tell()
                try:
                    r0 = next(sf)
                except StopIteration:
                    return cu_spills
                skip_qname = r0.query_name
                while True:
                    rec_voff = sf.tell()
                    try:
                        r = next(sf)
                    except StopIteration:
                        return cu_spills
                    if r.query_name != skip_qname:
                        # If `skip_qname` spanned all the way past `end_voff`, an earlier worker
                        # already processed `skip_qname` and the next worker starts at `r.query_name`.
                        if end_voff is not None and rec_voff > end_voff:
                            return cu_spills
                        first_rec = r
                        first_rec_voff = rec_voff
                        break
            else:
                rec_voff = sf.tell()
                try:
                    first_rec = next(sf)
                    first_rec_voff = rec_voff
                except StopIteration:
                    return cu_spills

            cur_qname = first_rec.query_name
            cur_group = [first_rec]
            passed_end = end_voff is not None and first_rec_voff >= end_voff

            while True:
                rec_voff = sf.tell()
                try:
                    r = next(sf)
                except StopIteration:
                    break

                qname = r.query_name
                if qname == cur_qname:
                    cur_group.append(r)
                    if end_voff is not None and rec_voff >= end_voff:
                        passed_end = True
                else:
                    _process_read_group(
                        cur_group,
                        strand_check,
                        min_qual,
                        min_dist,
                        indel_flank,
                        snp_by_chrom,
                        mut_pos,
                        muts_rds,
                        tracks,
                        mut_types,
                        mut_type_set,
                        base_set,
                        base_to_muttypes,
                        wr,
                        cU,
                        freq,
                        track_files,
                    )
                    if mut_pos and len(cU) >= CU_SPILL_LIMIT:
                        cu_spills.append(_flush_cu_spill(cU, spill_dir, worker_idx, spill_idx))
                        spill_idx += 1

                    if passed_end or (end_voff is not None and rec_voff > end_voff):
                        cur_group = None
                        break

                    cur_qname = qname
                    cur_group = [r]
                    if end_voff is not None and rec_voff >= end_voff:
                        passed_end = True

            if cur_group:
                _process_read_group(
                    cur_group,
                    strand_check,
                    min_qual,
                    min_dist,
                    indel_flank,
                    snp_by_chrom,
                    mut_pos,
                    muts_rds,
                    tracks,
                    mut_types,
                    mut_type_set,
                    base_set,
                    base_to_muttypes,
                    wr,
                    cU,
                    freq,
                    track_files,
                )

        if mut_pos and cU:
            cu_spills.append(_flush_cu_spill(cU, spill_dir, worker_idx, spill_idx))

        return cu_spills
    finally:
        if track_files is not None:
            for f in track_files:
                f.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Python implementation of TimeLapse / NR-seq mutation calling"
    )
    required_named = parser.add_argument_group("required named arguments")
    required_named.add_argument(
        "-b", "--bam", type=str, required=True, metavar="in_file.bam",
        help="Queryname-sorted BAM file to process",
    )
    parser.add_argument(
        "-o", "--output", type=str, default=None,
        help="Output _counts.csv or _counts.csv.gz path (default: <bam_prefix>_counts.csv)",
    )
    parser.add_argument(
        "--cUOutput", type=str, default=None,
        help="Output _cU.csv or _cU.csv.gz path when --mutPos is enabled (default: <bam_prefix>_cU.csv)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of parallel worker processes (default: 1)",
    )
    parser.add_argument(
        "--mutType", default="TC", type=str,
        help="Comma-separated mutation types to record (default: TC)",
    )
    parser.add_argument(
        "--reads", default="PE", type=str, choices=["PE", "SE"],
        help="Sequencing read layout: PE or SE (default: PE)",
    )
    parser.add_argument(
        "--minDist", default=-1, type=int, metavar="<int>",
        help="Base distance from read-end filter (default: -1)",
    )
    parser.add_argument(
        "--minQual", default=40, type=int, metavar="<int>",
        help="Base minimal quality filter (default: 40)",
    )
    parser.add_argument(
        "--indelFlank", default=3, type=int, metavar="<int>",
        help="Bases flanking any I/D CIGAR operation to mask from counts (default: 3)",
    )
    parser.add_argument(
        "--tracks", action="store_true",
        help="Generate read-name files for creating browser tracks",
    )
    parser.add_argument(
        "--mutsRDS", action="store_true",
        help="Track position-level mutation frequencies internally",
    )
    parser.add_argument(
        "--mutPos", action="store_true",
        help="Generate _cU file and per-read mutation position columns",
    )
    parser.add_argument("--SNPs", default=None, help="Path to snp.txt file")
    parser.add_argument(
        "--strandedness", default="F", type=str, choices=["F", "R"],
        help="Is first read forward or reverse orientation? F = forward is the default.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    mut_types = [m.strip() for m in args.mutType.split(",") if m.strip()]
    mut_type_set = set(mut_types)
    base_set = {m[0] for m in mut_types}
    base_to_muttypes = {}
    for mt in mut_types:
        base_to_muttypes.setdefault(mt[0], []).append(mt)

    input_name = args.bam[:-4] if args.bam.endswith(".bam") else args.bam.split(".bam")[0]
    counts_output = args.output if args.output else (input_name + "_counts.csv")
    cu_output = args.cUOutput if args.cUOutput else (input_name + "_cU.csv")
    is_gz = counts_output.endswith(".gz")

    header = ["qname", "nA", "nC", "nT", "nG", "rname", "FR", "sj"] + list(MUT_KEYS)
    if args.mutPos:
        header.extend(["gmutloc", "tp"])

    os.makedirs(os.path.dirname(os.path.abspath(counts_output)), exist_ok=True)
    if args.mutPos:
        os.makedirs(os.path.dirname(os.path.abspath(cu_output)), exist_ok=True)

    print("Start: " + str(datetime.datetime.now()))

    snp_by_chrom = load_snps(args.SNPs)
    effective_threads = 1 if (args.tracks or args.mutsRDS) else max(1, args.threads)
    chunks, is_empty = compute_bam_chunks(args.bam, effective_threads)

    if is_empty:
        opener = (
            (lambda p: gzip.open(p, "wt", newline="", compresslevel=3))
            if is_gz
            else (lambda p: open(p, "w", newline=""))
        )
        with opener(counts_output) as out_fh:
            csv.writer(out_fh).writerow(header)
        if args.mutPos:
            _merge_cu_spills([], cu_output)
        if args.tracks:
            for mt in mut_types:
                for i in range(6):
                    open(f"{input_name}_{mt}_{i}_reads.txt", "w").close()
        print("end: " + str(datetime.datetime.now()))
        return

    out_dir = os.path.dirname(os.path.abspath(counts_output))
    tmp_dir = tempfile.mkdtemp(prefix="mut_call_", dir=out_dir)
    try:
        cfg = {
            "bam_path": args.bam,
            "header": header,
            "is_gz": is_gz,
            "strand_check": args.strandedness == "F",
            "min_qual": args.minQual,
            "min_dist": args.minDist,
            "indel_flank": args.indelFlank,
            "snp_by_chrom": snp_by_chrom,
            "mut_pos": args.mutPos,
            "muts_rds": args.mutsRDS,
            "tracks": args.tracks,
            "mut_types": mut_types,
            "mut_type_set": mut_type_set,
            "base_set": base_set,
            "base_to_muttypes": base_to_muttypes,
            "spill_dir": tmp_dir,
            "track_prefix": input_name,
        }

        if len(chunks) == 1:
            _init_worker(cfg)
            all_cu_spills = _run_chunk_worker((0, chunks[0][0], chunks[0][1], counts_output))
        else:
            tasks = []
            chunk_paths = []
            ext = ".csv.gz" if is_gz else ".csv"
            for idx, (s_voff, e_voff) in enumerate(chunks):
                cpath = os.path.join(tmp_dir, f"chunk_{idx:04d}{ext}")
                chunk_paths.append(cpath)
                tasks.append((idx, s_voff, e_voff, cpath))

            all_cu_spills = []
            with mp.Pool(processes=len(chunks), initializer=_init_worker, initargs=(cfg,)) as pool:
                for worker_spills in pool.map(_run_chunk_worker, tasks):
                    all_cu_spills.extend(worker_spills)

            with open(counts_output, "wb") as dst_fh:
                for cpath in chunk_paths:
                    if os.path.exists(cpath):
                        with open(cpath, "rb") as src_fh:
                            shutil.copyfileobj(src_fh, dst_fh, length=16 * 1024 * 1024)
                        os.remove(cpath)

        print("end: " + str(datetime.datetime.now()))

        if args.mutPos:
            _merge_cu_spills(all_cu_spills, cu_output)
            print("cU: " + str(datetime.datetime.now()))
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
