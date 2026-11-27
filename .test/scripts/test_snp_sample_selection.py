#!/usr/bin/env python3
"""Regression tests for no-control SNP sample selection and depth ranking."""

import ast
import csv
import gzip
import importlib.util
import os
import random
import struct
import subprocess
import sys
import tempfile


BASE_CFG = {
    "call_snps_without_ctl": True,
    "snp_min_depth": 5,
    "snp_min_freq": 0.25,
    "snp_err_rate": 0.002,
    "snp_noctl_min_depth": 20,
    "snp_noctl_min_freq": 0.35,
    "snp_noctl_err_rate": 0.10,
    "snp_max_samples": -1,
    "snp_samples": [],
    "snp_sample_strategy": "random",
    "snp_random_seed": 42,
}


def _load_module(module_path, name):
    spec = importlib.util.spec_from_file_location(name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_common_helpers(common_smk_path):
    """Extract and compile SNP sample selection helpers from common.smk."""
    with open(common_smk_path, "r", encoding="utf-8") as fh:
        source = fh.read()
    tree = ast.parse(source, filename=common_smk_path)
    target_names = {"_parse_snp_samples", "resolve_snp_sample_selection"}
    kept_nodes = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in target_names
    ]
    assert len(kept_nodes) == 2, f"Expected 2 helper defs, found {len(kept_nodes)}"
    mod_ast = ast.Module(body=kept_nodes, type_ignores=[])
    namespace = {"random": random}
    exec(compile(mod_ast, common_smk_path, "exec"), namespace)
    return namespace["resolve_snp_sample_selection"]


def _write_synthetic_bam(path, n_records):
    """Write a minimal spec-valid BAM file containing `n_records` dummy records."""
    sam_hdr = b"@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n"
    ref_name = b"chr1\x00"
    # Minimal 36-byte BAM alignment block (read_name = "r\x00", 0 cigar, 0 seq)
    rec_payload = struct.pack(
        "<iiIIiiii",
        0,  # refID
        10,  # pos
        (2 << 0) | (20 << 8) | (4680 << 16),  # l_read_name=2, mapq=20, bin=4680
        0,  # flag_nc = 0
        0,  # l_seq = 0
        -1,  # next_refID
        -1,  # next_pos
        0,  # tlen
    ) + b"r\x00"
    rec_bytes = struct.pack("<I", len(rec_payload)) + rec_payload

    with gzip.open(path, "wb") as gz_fh:
        gz_fh.write(b"BAM\x01")
        gz_fh.write(struct.pack("<I", len(sam_hdr)) + sam_hdr)
        gz_fh.write(struct.pack("<I", 1))
        gz_fh.write(struct.pack("<I", len(ref_name)) + ref_name)
        gz_fh.write(struct.pack("<I", 1000))
        for _ in range(n_records):
            gz_fh.write(rec_bytes)


def _expect_value_error(fn, *args):
    try:
        fn(*args)
    except ValueError:
        return
    raise AssertionError(f"Expected ValueError for args={args}")


def run_all_tests(selector_path, common_smk_path, test_bam_dir):
    selector = _load_module(selector_path, "select_snp_samples")
    resolve_cfg = _load_common_helpers(common_smk_path)

    samps = ["S1", "S2", "S3", "S4", "S5"]

    # 1. Control samples present -> ignores no-ctl subsampling & uses standard cutoffs
    cfg_ctl = dict(
        BASE_CFG,
        snp_max_samples=1,
        snp_samples=["S3"],
        snp_sample_strategy="depth",
    )
    res_ctl = resolve_cfg(cfg_ctl, samps, ["S1", "S2"])
    assert res_ctl["use_depth_checkpoint"] is False
    assert res_ctl["selected_samps"] == ["S1", "S2"]
    assert res_ctl["min_depth"] == 5

    # 2. No controls, call_snps_without_ctl=False -> empty selection
    cfg_off = dict(BASE_CFG, call_snps_without_ctl=False, snp_max_samples=2)
    res_off = resolve_cfg(cfg_off, samps, [])
    assert res_off["use_depth_checkpoint"] is False
    assert res_off["selected_samps"] == []

    # 3. Default snp_max_samples=-1 -> uses all samples statically
    res_all = resolve_cfg(dict(BASE_CFG), samps, [])
    assert res_all["use_depth_checkpoint"] is False
    assert res_all["selected_samps"] == samps
    assert res_all["min_depth"] == 20

    # 4. snp_max_samples >= len(candidates) -> no checkpoint even when strategy="depth"
    cfg_over = dict(BASE_CFG, snp_max_samples=10, snp_sample_strategy="depth")
    res_over = resolve_cfg(cfg_over, samps, [])
    assert res_over["use_depth_checkpoint"] is False
    assert res_over["selected_samps"] == samps

    # 5. Explicit snp_samples (list and CSV string, deduplicated, order-preserving)
    cfg_exp = dict(BASE_CFG, snp_samples=["S4", "S2", "S4"])
    res_exp = resolve_cfg(cfg_exp, samps, [])
    assert res_exp["selected_samps"] == ["S4", "S2"]

    cfg_csv = dict(BASE_CFG, snp_samples="S5, S1, S3", snp_max_samples=2)
    res_csv = resolve_cfg(cfg_csv, samps, [])
    assert len(res_csv["selected_samps"]) == 2
    assert set(res_csv["selected_samps"]).issubset({"S5", "S1", "S3"})

    # 6. Seeded random determinism, dict-order set invariance, global RNG isolation
    random.seed(999)
    before_rand = random.random()
    random.seed(999)
    cfg_rand = dict(
        BASE_CFG, snp_max_samples=3, snp_sample_strategy="random", snp_random_seed=7
    )
    res_r1 = resolve_cfg(cfg_rand, samps, [])
    after_rand = random.random()
    assert before_rand == after_rand, "Global random state was mutated"

    res_r2 = resolve_cfg(cfg_rand, list(reversed(samps)), [])
    assert len(res_r1["selected_samps"]) == 3
    assert set(res_r1["selected_samps"]) == set(res_r2["selected_samps"])
    # Preserves input sample relative order
    assert res_r1["selected_samps"] == [
        s for s in samps if s in set(res_r1["selected_samps"])
    ]

    # 7. Strategy == "depth" with 1 <= max_samples < len(candidates) triggers checkpoint
    cfg_dep = dict(
        BASE_CFG,
        snp_samples=["S2", "S3", "S5"],
        snp_max_samples=2,
        snp_sample_strategy="depth",
    )
    res_dep = resolve_cfg(cfg_dep, samps, [])
    assert res_dep["use_depth_checkpoint"] is True
    assert res_dep["selected_samps"] == []
    assert res_dep["candidate_samps"] == ["S2", "S3", "S5"]

    # 8. Invalid config guardrails
    for bad_max in (0, -2, True, False, "zero"):
        _expect_value_error(
            resolve_cfg, dict(BASE_CFG, snp_max_samples=bad_max), samps, []
        )
    _expect_value_error(
        resolve_cfg, dict(BASE_CFG, snp_sample_strategy="most_reads"), samps, []
    )
    _expect_value_error(
        resolve_cfg, dict(BASE_CFG, snp_random_seed=True), samps, []
    )
    _expect_value_error(
        resolve_cfg, dict(BASE_CFG, snp_samples=["S1", "MISSING"]), samps, []
    )

    # 9. BAM record counter & depth selector on real test BAMs + synthetic BAM fixtures
    real_bams = [
        os.path.join(test_bam_dir, name)
        for name in ("WT_nos4U.bam", "WT_replicate_1.bam", "WT_replicate_2.bam")
    ]
    for bam_file in real_bams:
        py_cnt = selector.count_bam_records(bam_file, use_samtools=False)
        auto_cnt = selector.count_bam_records(bam_file, use_samtools=True)
        assert py_cnt > 0 and py_cnt == auto_cnt

    with tempfile.TemporaryDirectory(prefix="snp_sel_test_") as tmp_dir:
        fixture_specs = [
            ("LowCov", 2),
            ("TieSampleB", 7),
            ("TieSampleA", 7),
            ("HighCov.s", 12),
        ]
        fixture_bams = []
        for s_name, n_rec in fixture_specs:
            bam_p = os.path.join(tmp_dir, f"{s_name}.s.bam")
            _write_synthetic_bam(bam_p, n_rec)
            fixture_bams.append(bam_p)

        sel_top2, ranked = selector.rank_and_select_samples(
            fixture_bams, max_samples=2, threads=2, use_samtools=False
        )
        assert sel_top2 == ["HighCov.s", "TieSampleA"]
        assert ranked == [
            ("HighCov.s", 12, True),
            ("TieSampleA", 7, True),
            ("TieSampleB", 7, False),
            ("LowCov", 2, False),
        ]

        sel_capped_high, _ = selector.rank_and_select_samples(
            fixture_bams, max_samples=10, threads=2, use_samtools=False
        )
        assert len(sel_capped_high) == 4

        out_sel = os.path.join(tmp_dir, "snps", "selected_samples.txt")
        out_dep = os.path.join(tmp_dir, "snps", "sample_depths.tsv")
        proc = subprocess.run(
            [
                sys.executable,
                selector_path,
                "--bams",
                *fixture_bams,
                "--max-samples",
                "2",
                "--threads",
                "2",
                "--output-selected",
                out_sel,
                "--output-depths",
                out_dep,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, f"CLI failed: {proc.stderr}"
        with open(out_sel, "r", encoding="utf-8") as fh:
            cli_selected = [ln.strip() for ln in fh if ln.strip()]
        assert cli_selected == ["HighCov.s", "TieSampleA"]

        with open(out_dep, "r", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [r["sample"] for r in rows] == [
            "HighCov.s",
            "TieSampleA",
            "TieSampleB",
            "LowCov",
        ]
        assert [r["selected"] for r in rows] == ["True", "True", "False", "False"]

    # 10. Guard workflow Snakefile & .smk rules against PEP 701 f-string regressions
    #     and verify download_fastqs + skip_trimming helper paths
    _assert_no_fstrings_in_workflow_smk(common_smk_path)
    _verify_sra_fastq_helpers(common_smk_path)


def _find_fstrings_in_file(path):
    """Return 1-based line numbers containing f-string literals in `path`."""
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()

    hits = []
    idx = 0
    length = len(text)
    line_no = 1
    in_str = None

    while idx < length:
        ch = text[idx]
        if ch == "\n":
            line_no += 1
            idx += 1
            continue

        if in_str is not None:
            if ch == "\\":
                if idx + 1 < length and text[idx + 1] == "\n":
                    line_no += 1
                idx += 2
                continue
            if text.startswith(in_str, idx):
                idx += len(in_str)
                in_str = None
                continue
            idx += 1
            continue

        if ch == "#":
            while idx < length and text[idx] != "\n":
                idx += 1
            continue

        prev_is_ident = idx > 0 and (
            text[idx - 1].isalnum() or text[idx - 1] == "_"
        )
        if not prev_is_ident and ch in "fFrRbBuU":
            p_end = idx
            while p_end < length and text[p_end] in "fFrRbBuU":
                p_end += 1
            if p_end < length and text[p_end] in ("'", '"'):
                prefix = text[idx:p_end].lower()
                q_char = text[p_end]
                delim = (
                    q_char * 3
                    if text.startswith(q_char * 3, p_end)
                    else q_char
                )
                if "f" in prefix:
                    hits.append(line_no)
                in_str = delim
                idx = p_end + len(delim)
                continue

        if ch in ("'", '"'):
            delim = ch * 3 if text.startswith(ch * 3, idx) else ch
            in_str = delim
            idx += len(delim)
            continue

        idx += 1

    return hits


def _assert_no_fstrings_in_workflow_smk(common_smk_path):
    """Ensure Snakefile and workflow/rules/*.smk stay free of f-strings."""
    rules_dir = os.path.dirname(os.path.abspath(common_smk_path))
    workflow_dir = os.path.dirname(rules_dir)
    if not os.path.exists(os.path.join(workflow_dir, "Snakefile")):
        alt_workflow = os.path.abspath("../workflow")
        if os.path.exists(os.path.join(alt_workflow, "Snakefile")):
            workflow_dir = alt_workflow
            rules_dir = os.path.join(workflow_dir, "rules")

    smk_files = [os.path.join(workflow_dir, "Snakefile")]
    for entry in sorted(os.listdir(rules_dir)):
        if entry.endswith(".smk"):
            smk_files.append(os.path.join(rules_dir, entry))
    assert len(smk_files) >= 10, f"Expected workflow .smk files, got {smk_files}"

    for smk_file in smk_files:
        bad_lines = _find_fstrings_in_file(smk_file)
        assert not bad_lines, (
            f"Found f-string(s) in {smk_file} at line(s) {bad_lines}; "
            "use .format() or + in .smk files for Python 3.12 + Snakemake 7 safety."
        )

    # Sanity-check the lexer itself detects f-strings when present
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".smk", encoding="utf-8", delete=False
    ) as tmp_fh:
        tmp_fh.write('# f"ignored comment"\nx = "RF"\ny = f"bad_{x}"\n')
        probe_path = tmp_fh.name
    try:
        assert _find_fstrings_in_file(probe_path) == [3]
    finally:
        os.remove(probe_path)


def _verify_sra_fastq_helpers(common_smk_path):
    """Verify SRA download_fastqs path resolution with and without skip_trimming."""
    with open(common_smk_path, "r", encoding="utf-8") as fh:
        source = fh.read()
    tree = ast.parse(source, filename=common_smk_path)
    wanted = {
        "_list_sample_fastqs",
        "get_input_fastqs",
        "get_fastqc_read",
        "get_fastq_r1",
        "get_fastq_r2",
        "get_hisat2_reads",
    }
    fn_nodes = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    assert len(fn_nodes) == len(wanted)

    class _DummyLogger:
        def debug(self, *_args):
            pass

    class _WC:
        def __init__(self, sample, read="1"):
            self.sample = sample
            self.read = read

    cfg = {
        "download_fastqs": True,
        "PE": True,
        "skip_trimming": True,
        "samples": {},
    }
    ns = {
        "os": os,
        "config": cfg,
        "logger": _DummyLogger(),
        "is_gz": False,
        "READS": [1, 2],
        "VALID_FASTQ_EXTS": (".fastq", ".fq", ".fastq.gz", ".fq.gz"),
        "expand": lambda pat, **kw: [pat],
    }
    compiled = compile(
        ast.Module(body=fn_nodes, type_ignores=[]), common_smk_path, "exec"
    )
    exec(compiled, ns)

    wc1 = _WC("SRR24069126", "1")
    wc2 = _WC("SRR24069126", "2")
    expected_pe = [
        "results/download_fastq/SRR24069126_1.fastq",
        "results/download_fastq/SRR24069126_2.fastq",
    ]
    assert ns["get_input_fastqs"](wc1) == expected_pe
    assert ns["get_fastq_r1"](wc1) == expected_pe[0]
    assert ns["get_fastq_r2"](wc1) == expected_pe[1]
    assert ns["get_fastqc_read"](wc1) == expected_pe[0]
    assert ns["get_fastqc_read"](wc2) == expected_pe[1]
    assert ns["get_hisat2_reads"](wc1, READS=[1, 2]) == expected_pe

    cfg["PE"] = False
    expected_se = ["results/download_fastq/SRR24069126.fastq"]
    assert ns["get_input_fastqs"](wc1) == expected_se
    assert ns["get_fastq_r1"](wc1) == expected_se[0]
    assert ns["get_fastqc_read"](wc1) == expected_se[0]
    assert ns["get_hisat2_reads"](wc1, READS=[1]) == expected_se


def main():
    if len(sys.argv) != 5:
        sys.stderr.write(
            "Usage: test_snp_sample_selection.py <select_snp_samples.py> "
            "<common.smk> <bams_dir> <sentinel.ok>\n"
        )
        sys.exit(1)

    selector_path, common_smk_path, test_bam_dir, sentinel_out = sys.argv[1:5]
    run_all_tests(selector_path, common_smk_path, test_bam_dir)

    out_dir = os.path.dirname(os.path.abspath(sentinel_out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(sentinel_out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("ok\n")
    print("All SNP sample selection tests passed.")


if __name__ == "__main__":
    main()
