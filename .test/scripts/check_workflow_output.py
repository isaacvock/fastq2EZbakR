#!/usr/bin/env python3
"""Validate end-to-end fastq2EZbakR workflow outputs and clean run artifacts."""

import argparse
import csv
import glob
import gzip
import os
import shutil
import sys
import tempfile


FEATURE_COLUMNS = (
    ("genes", ("GF",), True),
    ("exons", ("XF",), True),
    ("transcripts", ("transcripts",), False),
    ("exonic_bins", ("exon_bin",), False),
    ("tec", ("TEC",), False),
    ("junctions", ("junction_start", "junction_end"), False),
    ("threeputr", ("threepUTR",), False),
    ("eej", ("ee_junction_id",), False),
    ("eij", ("ei_junction_id",), False),
)


def _strip_yaml_comment(line):
    """Remove trailing `# ...` comments outside single or double quotes."""
    in_quote = None
    escape = False
    for idx, ch in enumerate(line):
        if escape:
            escape = False
            continue
        if ch == "\\" and in_quote == '"':
            escape = True
            continue
        if ch in ("'", '"'):
            if in_quote is None:
                in_quote = ch
            elif in_quote == ch:
                in_quote = None
        elif ch == "#" and in_quote is None:
            return line[:idx]
    return line


def _parse_scalar(raw):
    val = raw.strip()
    if val in ("True", "true"):
        return True
    if val in ("False", "false"):
        return False
    if val in ("None", "null", "~"):
        return None
    if val.startswith("[") and val.endswith("]"):
        inner = val[1:-1].strip()
        if not inner:
            return []
        return [_parse_scalar(part) for part in inner.split(",") if part.strip()]
    if len(val) >= 2 and val[0] == val[-1] and val[0] in ("'", '"'):
        inner = val[1:-1]
        if val[0] == '"':
            inner = inner.replace('\\"', '"').replace("\\\\", "\\")
        return inner
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        return val


def load_simple_yaml(path):
    """Parse 1- and 2-level Snakemake YAML config files using stdlib only."""
    cfg = {}
    cur_section = None
    with open(path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            cleaned = _strip_yaml_comment(raw_line).rstrip()
            if not cleaned.strip():
                continue
            indent = len(cleaned) - len(cleaned.lstrip(" "))
            key_part, sep, val_part = cleaned.strip().partition(":")
            if not sep:
                continue
            key = key_part.strip()
            val_str = val_part.strip()
            if indent == 0:
                if val_str == "":
                    cfg[key] = {}
                    cur_section = key
                else:
                    cfg[key] = _parse_scalar(val_str)
                    cur_section = None
            elif cur_section is not None:
                cfg[cur_section][key] = _parse_scalar(val_str)
    return cfg


def _expected_columns(cfg):
    feats = cfg.get("features", {})
    feat_cols = []
    for feat_key, cols, default_val in FEATURE_COLUMNS:
        if feats.get(feat_key, default_val):
            feat_cols.extend(cols)

    mut_raw = str(cfg.get("mut_tracks", "TC"))
    mut_cols = [m.strip() for m in mut_raw.split(",") if m.strip()]
    base_cols = list(dict.fromkeys("n" + m[0] for m in mut_cols))
    return feat_cols, mut_cols, base_cols


def _check_gz_csv(path, expected_header, expected_samples, count_cols=(), min_n=1):
    assert os.path.isfile(path), f"Expected output file missing: {path}"
    seen_samples = set()
    row_count = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as gz_fh:
        reader = csv.reader(gz_fh)
        header = next(reader, None)
        assert header == expected_header, (
            f"Header mismatch in {path}: got {header}, expected {expected_header}"
        )
        col_idx = {name: idx for idx, name in enumerate(header)}
        for row in reader:
            assert len(row) == len(expected_header), f"Malformed row in {path}: {row}"
            row_count += 1
            seen_samples.add(row[col_idx["sample"]])
            n_val = int(row[col_idx["n"]])
            assert n_val >= min_n, (
                f"Expected n >= {min_n} in {path}, got {n_val}"
            )
            for c_name in count_cols:
                c_val = int(row[col_idx[c_name]])
                assert c_val >= 0, f"Negative count in {path} col {c_name}: {c_val}"

    assert row_count > 0, f"Expected >=1 data row in {path}, got 0"
    missing_samps = set(expected_samples) - seen_samples
    assert not missing_samps, (
        f"Samples {sorted(missing_samps)} missing from {path} (saw {seen_samples})"
    )


def _check_parquet_file(path, required_col_names):
    assert os.path.isfile(path), f"Missing Parquet partition: {path}"
    size = os.path.getsize(path)
    assert size > 12, f"Parquet file too small ({size} bytes): {path}"
    with open(path, "rb") as fh:
        data = fh.read()
    assert data[:4] == b"PAR1" and data[-4:] == b"PAR1", (
        f"Invalid Parquet magic bytes in {path}"
    )
    for col_name in required_col_names:
        assert col_name.encode("utf-8") in data, (
            f"Column '{col_name}' not found in Parquet schema metadata of {path}"
        )


def validate_workflow_outputs(config_path, workdir="."):
    """Assert that all expected workflow outputs in `workdir` match `config_path`."""
    cfg = load_simple_yaml(config_path)
    samples_map = cfg.get("samples", {})
    download_fastqs = bool(cfg.get("download_fastqs", False))
    if download_fastqs:
        samp_names = [str(s).strip() for s in cfg.get("sra_accessions", [])]
    else:
        samp_names = list(samples_map.keys())
    assert samp_names, f"No sample names found in {config_path}"

    bam2bakr = (not download_fastqs) and all(
        str(v).endswith(".bam") for v in samples_map.values()
    )
    aligner = str(cfg.get("aligner", "star")).strip().lower().replace("-", "")
    final_out = cfg.get(
        "final_output", {"cB": True, "cUP": False, "arrow": False}
    )
    make_cb = bool(final_out.get("cB", False))
    make_cup = bool(final_out.get("cUP", False))
    make_arrow = bool(final_out.get("arrow", False))
    make_tracks = bool(cfg.get("make_tracks", True))
    mutpos = bool(cfg.get("mutpos", False))
    run_rsem = bool(cfg.get("run_rsem", True))

    feat_cols, mut_cols, base_cols = _expected_columns(cfg)

    def _wpath(*parts):
        return os.path.join(workdir, *parts)

    cb_path = _wpath("results", "cB", "cB.csv.gz")
    if make_cb:
        cb_header = ["sample", "rname", "sj", *feat_cols, *mut_cols, *base_cols, "n"]
        _check_gz_csv(
            cb_path, cb_header, samp_names, count_cols=(*mut_cols, *base_cols)
        )
    else:
        assert not os.path.exists(cb_path), (
            f"Unexpected cB output when cB=False: {cb_path}"
        )

    cup_path = _wpath("results", "cUP", "cUP.csv.gz")
    if make_cup:
        cup_header = ["sample", "rname", *feat_cols, *mut_cols, *base_cols, "n"]
        _check_gz_csv(cup_path, cup_header, samp_names, count_cols=mut_cols)
    else:
        assert not os.path.exists(cup_path), (
            f"Unexpected cUP output when cUP=False: {cup_path}"
        )

    arrow_dir = _wpath("results", "arrow_dataset")
    if make_arrow:
        req_cols = ["sample", "rname", "sj", *feat_cols, *mut_cols, *base_cols, "n"]
        for s_name in samp_names:
            pq_path = os.path.join(arrow_dir, f"sample={s_name}", "part-0.parquet")
            _check_parquet_file(pq_path, req_cols)
    else:
        assert not os.path.exists(arrow_dir), (
            f"Unexpected arrow_dataset when arrow=False: {arrow_dir}"
        )

    if mutpos:
        mutpos_hdr = ["sample", "rname", "gloc", "tp", "trials", "n"]
        for fname in ("mutpos.csv.gz", "mutpos_filtered.csv.gz"):
            _check_gz_csv(
                _wpath("results", "cB", fname),
                mutpos_hdr,
                samp_names,
                count_cols=("trials",),
                min_n=0,
            )

    tracks_dir = _wpath("results", "tracks")
    if make_tracks:
        for s_name in samp_names:
            for mut in mut_cols:
                for level in range(6):
                    for strand in ("pos", "min"):
                        tdf = os.path.join(
                            tracks_dir, f"{s_name}.{mut}.{level}.{strand}.tdf"
                        )
                        assert os.path.isfile(tdf) and os.path.getsize(tdf) > 0, (
                            f"Missing or empty TDF track: {tdf}"
                        )
    else:
        assert not os.path.exists(tracks_dir), (
            f"Unexpected tracks directory when make_tracks=False: {tracks_dir}"
        )

    if not bam2bakr:
        mqc = _wpath("results", "multiqc", "multiqc_report.html")
        assert os.path.isfile(mqc) and os.path.getsize(mqc) > 0, (
            f"Missing MultiQC report: {mqc}"
        )
        if aligner == "star" and run_rsem:
            for s_name in samp_names:
                iso = _wpath("results", "rsem", f"{s_name}.isoforms.results")
                assert os.path.isfile(iso) and os.path.getsize(iso) > 0, (
                    f"Missing RSEM isoforms output: {iso}"
                )

    if cfg.get("features", {}).get("exonic_bins", False):
        flat_gtf = _wpath(str(cfg.get("flat_annotation", "data/annotation/flat.gtf")))
        assert os.path.isfile(flat_gtf) and os.path.getsize(flat_gtf) > 0, (
            f"Missing flattened GTF: {flat_gtf}"
        )
        with open(flat_gtf, "r", encoding="utf-8") as fh:
            flat_head = fh.read(4096)
        assert "exonic_part" in flat_head and "exon_id" in flat_head

    snp_txt = _wpath("results", "snps", "snp.txt")
    snp_tsv = _wpath("results", "snps", "snp_counts.tsv")
    assert os.path.isfile(snp_txt), f"Missing snp.txt: {snp_txt}"
    assert os.path.isfile(snp_tsv), f"Missing snp_counts.tsv: {snp_tsv}"

    ctl_samples = cfg.get("control_samples", [])
    noctl = bool(cfg.get("call_snps_without_ctl", False))
    max_samps = int(cfg.get("snp_max_samples", -1))
    strategy = str(cfg.get("snp_sample_strategy", "random")).strip().lower()
    explicit_pool = cfg.get("snp_samples", []) or samp_names

    if not ctl_samples and noctl and 1 <= max_samps < len(explicit_pool):
        if strategy == "depth":
            sel_txt = _wpath("results", "snps", "selected_samples.txt")
            dep_tsv = _wpath("results", "snps", "sample_depths.tsv")
            assert os.path.isfile(sel_txt), f"Missing checkpoint output: {sel_txt}"
            assert os.path.isfile(dep_tsv), f"Missing sample depths TSV: {dep_tsv}"
            with open(sel_txt, "r", encoding="utf-8") as fh:
                chosen = [ln.strip() for ln in fh if ln.strip()]
            assert len(chosen) == max_samps, (
                f"Expected {max_samps} depth-selected samples, got {chosen}"
            )
            assert set(chosen).issubset(set(explicit_pool))

        if not make_tracks and not cfg.get("features", {}).get("threeputr", False):
            sort_log_dir = _wpath("logs", "sort_coord_bam")
            if os.path.isdir(sort_log_dir):
                sort_logs = [
                    fn for fn in os.listdir(sort_log_dir) if fn.endswith(".log")
                ]
                assert len(sort_logs) == max_samps, (
                    f"Expected sort_coord_bam to run on {max_samps} samples, "
                    f"found logs: {sort_logs}"
                )

    print(f"Validated all workflow outputs for {config_path}.", flush=True)
    return cfg


def clean_run_artifacts(cfg, workdir=".", clean_indices=False):
    """Remove per-run outputs so subsequent test configs execute from scratch."""
    for rel_dir in ("results", "logs", "junction_annotation"):
        target = os.path.join(workdir, rel_dir)
        if os.path.exists(target):
            shutil.rmtree(target)

    for pattern in (
        os.path.join(workdir, "data", "annotation", "flat*.gtf"),
        os.path.join(workdir, "data", "genome", "*.fai"),
    ):
        for match in glob.glob(pattern):
            if os.path.isfile(match):
                os.remove(match)

    if clean_indices:
        idx_dir = str(cfg.get("indices", "")).strip()
        for rel_idx in filter(None, (idx_dir, "rsem_index")):
            target = os.path.join(workdir, rel_idx)
            if os.path.exists(target):
                shutil.rmtree(target)

    print(
        f"Cleaned run artifacts in {workdir} (clean_indices={clean_indices}).",
        flush=True,
    )


def _run_self_test(configs_dir):
    """Unit-test YAML parsing on all repo test configs and mock validation."""
    expected_configs = (
        "transcript_config.yaml",
        "star_rsem_config.yaml",
        "hisat3n_config.yaml",
        "bam2bakr_config.yaml",
        "hisat2_mutpos_config.yaml",
        "threeputr_config.yaml",
    )
    for name in expected_configs:
        cpath = os.path.join(configs_dir, name)
        assert os.path.isfile(cpath), f"Missing test config: {cpath}"
        parsed = load_simple_yaml(cpath)
        assert "samples" in parsed and "features" in parsed, f"Bad parse: {name}"

    # Simulate HISAT-3N depth-subsampled output tree and verify validator + cleaner
    h3n_cfg_path = os.path.join(configs_dir, "hisat3n_config.yaml")
    with tempfile.TemporaryDirectory(prefix="wf_val_test_") as tmp_dir:
        os.makedirs(os.path.join(tmp_dir, "results", "cB"), exist_ok=True)
        os.makedirs(os.path.join(tmp_dir, "results", "multiqc"), exist_ok=True)
        os.makedirs(os.path.join(tmp_dir, "results", "snps"), exist_ok=True)
        os.makedirs(os.path.join(tmp_dir, "logs", "sort_coord_bam"), exist_ok=True)
        os.makedirs(os.path.join(tmp_dir, "data", "annotation"), exist_ok=True)
        os.makedirs(os.path.join(tmp_dir, "hisat3n_index"), exist_ok=True)

        cb_cols = [
            "sample",
            "rname",
            "sj",
            "GF",
            "XF",
            "exon_bin",
            "TEC",
            "junction_start",
            "junction_end",
            "TC",
            "GA",
            "nT",
            "nG",
            "n",
        ]
        with gzip.open(
            os.path.join(tmp_dir, "results", "cB", "cB.csv.gz"),
            "wt",
            encoding="utf-8",
            newline="",
        ) as fh:
            wr = csv.writer(fh)
            wr.writerow(cb_cols)
            for s_id in ("WT_1", "WT_2", "WT_3"):
                wr.writerow(
                    [
                        s_id,
                        "chr21",
                        "F",
                        "G1",
                        "G1",
                        "EG1001",
                        "T1",
                        "10",
                        "20",
                        1,
                        0,
                        15,
                        12,
                        4,
                    ]
                )

        for s_id in ("WT_1", "WT_2", "WT_3"):
            pq_dir = os.path.join(
                tmp_dir, "results", "arrow_dataset", f"sample={s_id}"
            )
            os.makedirs(pq_dir, exist_ok=True)
            meta = ",".join(cb_cols).encode("utf-8")
            with open(os.path.join(pq_dir, "part-0.parquet"), "wb") as pq_fh:
                pq_fh.write(b"PAR1" + meta + b"PAR1")

        with open(
            os.path.join(tmp_dir, "results", "multiqc", "multiqc_report.html"),
            "w",
            encoding="utf-8",
        ) as fh:
            fh.write("<html>ok</html>\n")

        with open(
            os.path.join(tmp_dir, "data", "annotation", "flat.gtf"),
            "w",
            encoding="utf-8",
        ) as fh:
            fh.write('chr21\tsrc\texonic_part\t1\t10\t.\t+\t.\texon_id "E1";\n')

        for fname in ("snp.txt", "snp_counts.tsv"):
            with open(
                os.path.join(tmp_dir, "results", "snps", fname),
                "w",
                encoding="utf-8",
            ) as fh:
                fh.write("\n")

        with open(
            os.path.join(tmp_dir, "results", "snps", "selected_samples.txt"),
            "w",
            encoding="utf-8",
        ) as fh:
            fh.write("WT_1\nWT_2\n")
        with open(
            os.path.join(tmp_dir, "results", "snps", "sample_depths.tsv"),
            "w",
            encoding="utf-8",
        ) as fh:
            fh.write("sample\tmapped_records\tselected\n")

        for s_id in ("WT_1", "WT_2"):
            with open(
                os.path.join(tmp_dir, "logs", "sort_coord_bam", f"{s_id}.log"),
                "w",
                encoding="utf-8",
            ) as fh:
                fh.write("sorted\n")

        cfg = validate_workflow_outputs(h3n_cfg_path, workdir=tmp_dir)
        clean_run_artifacts(cfg, workdir=tmp_dir, clean_indices=True)
        assert not os.path.exists(os.path.join(tmp_dir, "results"))
        assert not os.path.exists(os.path.join(tmp_dir, "hisat3n_index"))

    print("check_workflow_output.py self-tests passed.", flush=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Validate fastq2EZbakR test run outputs and optionally clean up."
    )
    parser.add_argument("--config", default="", help="Path to workflow config YAML.")
    parser.add_argument("--workdir", default=".", help="Test working directory.")
    parser.add_argument(
        "--clean",
        action="store_true",
        default=False,
        help="Remove results/logs/flat GTF after validation succeeds.",
    )
    parser.add_argument(
        "--clean-indices",
        action="store_true",
        default=False,
        help="Also remove aligner and RSEM index directories when --clean is set.",
    )
    parser.add_argument(
        "--output-sentinel",
        default="",
        help="Optional sentinel file to touch upon success.",
    )
    parser.add_argument(
        "--self-test",
        default="",
        help="Run built-in unit tests against the given .test/configs directory.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.self_test:
        _run_self_test(args.self_test)
    else:
        if not args.config:
            raise ValueError("--config is required unless --self-test is specified.")
        cfg = validate_workflow_outputs(args.config, workdir=args.workdir)
        if args.clean:
            clean_run_artifacts(
                cfg, workdir=args.workdir, clean_indices=args.clean_indices
            )

    if args.output_sentinel:
        out_dir = os.path.dirname(os.path.abspath(args.output_sentinel))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.output_sentinel, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("ok\n")


if __name__ == "__main__":
    main()
