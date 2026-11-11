import glob
import os


### Figure out if bam files or fastq files are provided
if config.get("download_fastqs", False):
    config["bam2bakr"] = False
else:
    _sample_vals = list(config.get("samples", {}).values())
    config["bam2bakr"] = bool(_sample_vals) and all(
        str(value).endswith(".bam") for value in _sample_vals
    )


### If assigning reads to junctions, need to keep jI and jM tags
if config["features"].get("junctions", False):
    config["remove_tags"] = False


### GENERAL HELPER FUNCTIONS/VARIABLES USED IN ALL CASES

if config["download_fastqs"]:
    SAMP_NAMES = list(config["sra_accessions"])
else:
    SAMP_NAMES = list(config["samples"].keys())

NUM_SAMPS = len(SAMP_NAMES)

# Directory containing index; used in case of certain aligners
INDEX_DIR = config["indices"]

# Strip trailing slash if provided
if str(config["indices"]).endswith("/"):
    INDEX_PATH = str(config["indices"])[:-1]
else:
    INDEX_PATH = str(config["indices"])

# Determine how many fastqs to look for
if config["PE"]:
    READS = [1, 2]
    READ_NAMES = ["r1", "r2"]
else:
    READS = [1]
    READ_NAMES = ["r1"]

VALID_FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def _list_sample_fastqs(fastq_path):
    if not os.path.isdir(fastq_path):
        raise ValueError(f"FASTQ directory does not exist: {fastq_path}")
    fastq_files = sorted(
        os.path.join(fastq_path, entry)
        for entry in os.listdir(fastq_path)
        if not entry.startswith(".") and entry.endswith(VALID_FASTQ_EXTS)
    )
    if not fastq_files:
        raise ValueError(
            f"No FASTQ files ({', '.join(VALID_FASTQ_EXTS)}) found in {fastq_path}"
        )
    return fastq_files


# Define the function to get input FASTQ files
# Suggested in issue #23: https://github.com/isaacvock/fastq2EZbakR/issues/23
def get_input_fastqs(wildcards):
    logger.debug(f"Resolving input for sample: {wildcards.sample}")

    if config.get("download_fastqs", False):
        sra_reads = ["_1.fastq", "_2.fastq"] if config.get("PE", False) else [".fastq"]
        paths = sorted(
            f"results/download_fastq/{wildcards.sample}{read}" for read in sra_reads
        )
        logger.debug(f"SRA FASTQ paths for {wildcards.sample}: {paths}")
        return paths

    fastq_path = config["samples"].get(wildcards.sample, None)
    if fastq_path is None:
        raise ValueError(f"No path found for sample {wildcards.sample}")

    fastq_files = _list_sample_fastqs(fastq_path)
    logger.debug(f"Found FASTQ files for {wildcards.sample}: {fastq_files}")
    return fastq_files


# Check if fastq files are gzipped
fastq_paths = config.get("samples", {})
is_gz = False

if not config["download_fastqs"] and not config["bam2bakr"]:
    for p in fastq_paths.values():
        fastqs = [
            path
            for path in set(glob.glob(f"{p}/*.fastq*") + glob.glob(f"{p}/*.fq*"))
            if path.endswith(VALID_FASTQ_EXTS)
        ]
        if any(path.endswith((".fastq.gz", ".fq.gz")) for path in fastqs):
            is_gz = True


# Get mutation types to track
MutTypes = config["mut_tracks"]
Mutation_Types = [m.strip() for m in MutTypes.split(",") if m.strip()]
Nucleotide_Types = [f"n{muttype[0]}" for muttype in Mutation_Types]


### FASTQC HELPERS


def get_fastqc_read(wildcards):
    if config["skip_trimming"]:
        if is_gz:
            return expand(
                "results/unzipped/{SID}.{READ}.fastq",
                SID=wildcards.sample,
                READ=wildcards.read,
            )
        fastq_path = config["samples"][wildcards.sample]
        fastq_files = _list_sample_fastqs(fastq_path)
        read_id = int(wildcards.read) - 1
        return fastq_files[read_id]

    return expand(
        "results/trimmed/{SID}.{READ}.fastq",
        SID=wildcards.sample,
        READ=wildcards.read,
    )


### STAR / HISAT-3N READ HELPERS


def get_fastq_r1(wildcards):
    if config["skip_trimming"]:
        if is_gz:
            return expand("results/unzipped/{SID}.1.fastq", SID=wildcards.sample)
        fastq_path = config["samples"][wildcards.sample]
        return _list_sample_fastqs(fastq_path)[0]

    return expand("results/trimmed/{SID}.1.fastq", SID=wildcards.sample)


def get_fastq_r2(wildcards):
    if config["skip_trimming"]:
        if is_gz:
            return expand("results/unzipped/{SID}.2.fastq", SID=wildcards.sample)
        fastq_path = config["samples"][wildcards.sample]
        return _list_sample_fastqs(fastq_path)[1]

    return expand("results/trimmed/{SID}.2.fastq", SID=wildcards.sample)


# Determine whether or not indexing needs to be done
star_indices = glob.glob(f"{INDEX_PATH}/genomeParameters.txt")

if config["aligner"] == "star":
    make_index = len(star_indices) == 0


### HISAT2 HELPERS


def get_hisat2_reads(wildcards, READS=READS):
    if config["skip_trimming"]:
        if is_gz:
            return expand(
                "results/unzipped/{SID}.{read}.fastq", SID=wildcards.sample, read=READS
            )
        fastq_path = config["samples"][wildcards.sample]
        return _list_sample_fastqs(fastq_path)

    return expand(
        "results/trimmed/{SID}.{read}.fastq", SID=wildcards.sample, read=READS
    )


hisat2_indices = glob.glob(f"{INDEX_PATH}/*.ht2*")

if config["aligner"] == "hisat2":
    make_index = len(hisat2_indices) == 0


### HISAT-3N HELPERS

hisat3n_indices = glob.glob(f"{INDEX_PATH}/*.3n.*.ht2*")

if config["aligner"] == "hisat3n":
    make_index = len(hisat3n_indices) == 0

if str(config.get("hisat3n_base_change", "")).strip():
    HISAT3N_BASE_CHANGE = str(config["hisat3n_base_change"]).strip()
else:
    if not Mutation_Types or len(Mutation_Types[0]) != 2:
        raise ValueError(
            f"Cannot derive hisat3n_base_change from mut_tracks='{config.get('mut_tracks')}'. "
            "Specify hisat3n_base_change explicitly (e.g. 'T,C')."
        )
    HISAT3N_BASE_CHANGE = f"{Mutation_Types[0][0]},{Mutation_Types[0][1]}"
    if config["aligner"] == "hisat3n" and len(Mutation_Types) > 1:
        logger.warning(
            f"HISAT-3N supports a single --base-change per index/alignment run; "
            f"using '{HISAT3N_BASE_CHANGE}' from first entry of mut_tracks='{config['mut_tracks']}'."
        )

if config.get("hisat3n_index_base"):
    HISAT3N_BASE = f"{INDEX_PATH}/{config['hisat3n_index_base']}"
elif hisat3n_indices:
    _detected_h3n_prefix = os.path.basename(sorted(hisat3n_indices)[0]).split(".3n.")[0]
    HISAT3N_BASE = f"{INDEX_PATH}/{_detected_h3n_prefix}"
else:
    HISAT3N_BASE = f"{INDEX_PATH}/hisat3n_index"

HISAT3N_REPEAT_FLAG = (
    "--repeat-index" if config.get("hisat3n_repeat_index", False) else ""
)

_h3n_user_tokens = str(config.get("hisat3n_align_params", "")).split()
_h3n_auto_flags = []
if "--rna-strandness" not in _h3n_user_tokens:
    if config["strandedness"] == "reverse":
        _h3n_auto_flags.append(
            "--rna-strandness RF" if config["PE"] else "--rna-strandness R"
        )
    elif config["strandedness"] in ("yes", "forward"):
        _h3n_auto_flags.append(
            "--rna-strandness FR" if config["PE"] else "--rna-strandness F"
        )
if (
    config.get("hisat3n_directional", True)
    and "--directional-mapping" not in _h3n_user_tokens
    and "--directional-mapping-reverse" not in _h3n_user_tokens
):
    if config["strandedness"] == "reverse":
        _h3n_auto_flags.append("--directional-mapping-reverse")
    elif config["strandedness"] in ("yes", "forward"):
        _h3n_auto_flags.append("--directional-mapping")

HISAT3N_ALIGN_EXTRA = " ".join(
    part
    for part in [
        " ".join(_h3n_auto_flags),
        str(config.get("hisat3n_align_params", "")).strip(),
    ]
    if part
)


### Extra parameters passed to STAR alignment

if config.get("quantify_premRNA", False):
    AandQ_ANNOTATION = "results/modify_annotation/modified_annotation.gtf"
else:
    AandQ_ANNOTATION = config["annotation"]

STAR_PARAMS = str(config.get("star_align_params", ""))
args = STAR_PARAMS.split()

if config["features"].get("junctions", False):
    tags = list(config["star_sam_tags"])
    if "jI" not in tags:
        tags.append("jI")
    if "jM" not in tags:
        tags.append("jM")
    sam_attributes = set(tags)
else:
    sam_attributes = set(config["star_sam_tags"])


def _find_flag_span(arg_list, flag):
    idx = arg_list.index(flag)
    next_flag_idx = len(arg_list)
    for i in range(idx + 1, len(arg_list)):
        if arg_list[i].startswith("--"):
            next_flag_idx = i
            break
    return idx, next_flag_idx


# Process --outSAMattributes (sorted for deterministic command line)
if "--outSAMattributes" in args:
    idx, next_flag_idx = _find_flag_span(args, "--outSAMattributes")
    existing_attributes = set(args[idx + 1 : next_flag_idx])
    args[idx + 1 : next_flag_idx] = sorted(existing_attributes.union(sam_attributes))
else:
    args.extend(["--outSAMattributes"] + sorted(sam_attributes))

# Process --quantMode
quant_mode = {"TranscriptomeSAM", "GeneCounts"}
if "--quantMode" in args:
    idx, next_flag_idx = _find_flag_span(args, "--quantMode")
    existing_modes = set(args[idx + 1 : next_flag_idx])
    args[idx + 1 : next_flag_idx] = sorted(existing_modes.union(quant_mode))
else:
    args.extend(["--quantMode"] + sorted(quant_mode))

# Force --outSAMtype to BAM Unsorted (skipping STAR coordinate sort before sort_filter)
if "--outSAMtype" in args:
    idx, next_flag_idx = _find_flag_span(args, "--outSAMtype")
    args[idx + 1 : next_flag_idx] = ["BAM", "Unsorted"]
else:
    args.extend(["--outSAMtype", "BAM", "Unsorted"])

# Force --sjdbGTFfile to AandQ_ANNOTATION without clobbering adjacent flags
if "--sjdbGTFfile" in args:
    idx, next_flag_idx = _find_flag_span(args, "--sjdbGTFfile")
    args[idx + 1 : next_flag_idx] = [str(AandQ_ANNOTATION)]
else:
    args.extend(["--sjdbGTFfile", str(AandQ_ANNOTATION)])

STAR_EXTRA = " ".join(args)


### RSEM STRANDEDNESS HELPER

_rsem_user_tokens = str(config.get("rsem_quant_params", "")).split()
if "--strandedness" in _rsem_user_tokens or "--forward-prob" in _rsem_user_tokens:
    RSEM_QUANT_EXTRA = str(config.get("rsem_quant_params", ""))
else:
    if config["strandedness"] == "reverse":
        _rsem_strand_flag = "--strandedness reverse"
    elif config["strandedness"] in ("yes", "forward"):
        _rsem_strand_flag = "--strandedness forward"
    else:
        _rsem_strand_flag = "--strandedness none"
    RSEM_QUANT_EXTRA = " ".join(
        part
        for part in [
            _rsem_strand_flag,
            str(config.get("rsem_quant_params", "")).strip(),
        ]
        if part
    )


### HISAT2 STRANDEDNESS & INDEX BASE

if config["strandedness"] == "no":
    HISAT2_STRANDEDNESS = ""
elif config["strandedness"] == "reverse":
    HISAT2_STRANDEDNESS = (
        "--rna-strandness RF" if config["PE"] else "--rna-strandness R"
    )
else:
    HISAT2_STRANDEDNESS = (
        "--rna-strandness FR" if config["PE"] else "--rna-strandness F"
    )

if config["hisat2_index_base"]:
    HISAT2_BASE = "{}/{}".format(INDEX_PATH, config["hisat2_index_base"])
else:
    HISAT2_BASE = "{}/{}".format(INDEX_PATH, "hisat2_index")


### BAM2EZBAKR & MUTATIONAL ANOMALY (SNP) PARAMETERS

CTL_NAMES = list(config["control_samples"])
s4U_SAMPS = list(set(SAMP_NAMES) - set(CTL_NAMES))
nctl = len(CTL_NAMES)
nsamps = len(SAMP_NAMES)

if nctl > 0:
    SNP_SAMP_NAMES = CTL_NAMES
    SNP_MIN_DEPTH = config["snp_min_depth"]
    SNP_MIN_FREQ = config["snp_min_freq"]
    SNP_ERR_RATE = config["snp_err_rate"]
elif config.get("call_snps_without_ctl", False):
    SNP_SAMP_NAMES = SAMP_NAMES
    SNP_MIN_DEPTH = config["snp_noctl_min_depth"]
    SNP_MIN_FREQ = config["snp_noctl_min_freq"]
    SNP_ERR_RATE = config["snp_noctl_err_rate"]
else:
    SNP_SAMP_NAMES = []
    SNP_MIN_DEPTH = config["snp_min_depth"]
    SNP_MIN_FREQ = config["snp_min_freq"]
    SNP_ERR_RATE = config["snp_err_rate"]

SNP_BAMS = expand("results/coord_bam/{s}.bam", s=SNP_SAMP_NAMES)
SNP_BAIS = expand("results/coord_bam/{s}.bam.bai", s=SNP_SAMP_NAMES)


def get_index_name():
    genome = config["genome"]
    return str(genome) + ".fai"


FORMAT = "PE" if config["PE"] else "SE"
NORMALIZE = config["normalize"]


def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


### TEC MODE HELPER

_raw_tec_mode = str(config.get("tec_mode", "auto")).strip().lower()
if _raw_tec_mode not in ("auto", "transcriptome", "genome"):
    raise ValueError(
        f"Unsupported tec_mode '{config.get('tec_mode')}'. Expected 'auto', 'transcriptome', or 'genome'."
    )
if _raw_tec_mode == "auto":
    TEC_MODE = (
        "transcriptome"
        if (config["aligner"] == "star" and not config["bam2bakr"])
        else "genome"
    )
elif _raw_tec_mode == "transcriptome":
    if config["aligner"] != "star" or config["bam2bakr"]:
        raise ValueError(
            "tec_mode='transcriptome' requires aligner='star' with FASTQ input (bam2bakr=False). "
            "Use tec_mode='auto' or tec_mode='genome'."
        )
    TEC_MODE = "transcriptome"
else:
    TEC_MODE = "genome"


### TRACK NORMALIZATION HELPERS

if config.get("use_exons_only", True):
    NORMALIZATION_INPUT = expand(
        "results/featurecounts_exons/{sample}.featureCounts", sample=SAMP_NAMES
    )
    NORMALIZATION_EXTRA = ""
else:
    NORMALIZATION_INPUT = expand(
        "results/featurecounts_genes/{sample}.featureCounts", sample=SAMP_NAMES
    )
    NORMALIZATION_EXTRA = "--use_genes"


### 3'-UTR HELPERS

if config.get("strandedness") == "reverse":
    FC_3UTR_PARAMS = "-R CORE -t 3UTR -O --read2pos 5 -g utr_id"
else:
    FC_3UTR_PARAMS = "-R CORE -t 3UTR -O --read2pos 3 -g utr_id"

FC_3UTR_CORE_EXT = (
    "_informative.bam.featureCounts"
    if config.get("PE", False)
    else ".s.bam.featureCounts"
)

if config.get("call_threeputrs", False):
    THREEPUTR_ANNOTATION = "annotations/threepUTR_annotation.gtf"
else:
    THREEPUTR_ANNOTATION = "annotations/annotated_threepUTR_annotation.gtf"


def fetch_informative_read(wildcards):
    if config.get("PE", False):
        return expand(
            "results/informative_read/{sample}_informative.bam",
            sample=wildcards.sample,
        )
    return expand("results/sf_reads/{sample}.s.bam", sample=wildcards.sample)


### FEATURECOUNTS & MERGE HELPERS


def get_merge_input(wildcards):
    merge_input = [f"results/counts/{wildcards.sample}_counts.csv.gz"]

    if config.get("features").get("genes", True):
        merge_input.append(
            f"results/featurecounts_genes/{wildcards.sample}.s.bam.featureCounts"
        )

    if config.get("features").get("exons", True):
        merge_input.append(
            f"results/featurecounts_exons/{wildcards.sample}.s.bam.featureCounts"
        )

    if config.get("features").get("transcripts", False):
        merge_input.append(
            f"results/featurecounts_transcripts/{wildcards.sample}.s.bam.featureCounts"
        )

    if config.get("features").get("exonic_bins", False):
        merge_input.append(
            f"results/featurecounts_exonbins/{wildcards.sample}.s.bam.featureCounts"
        )

    if config.get("features").get("tec", False):
        merge_input.append(f"results/read_to_transcripts/{wildcards.sample}.csv")

    if config.get("features").get("junctions", False):
        merge_input.append(f"results/read_to_junctions/{wildcards.sample}.csv.gz")

    if config.get("features").get("threeputr", False):
        merge_input.append(
            f"results/featurecounts_3utr/{wildcards.sample}{FC_3UTR_CORE_EXT}"
        )

    if config.get("features").get("eej", False):
        merge_input.append(
            f"results/featurecounts_eej/{wildcards.sample}.s.bam.featureCounts"
        )

    if config.get("features").get("eij", False):
        merge_input.append(
            f"results/featurecounts_eij/{wildcards.sample}.s.bam.featureCounts"
        )

    return merge_input


if config["strandedness"] == "reverse":
    FC_STRAND = 2
elif config["strandedness"] == "yes" or config["strandedness"] == "forward":
    FC_STRAND = 1
else:
    FC_STRAND = 0

if config["PE"]:
    FC_GENES_PARAMS = " -R CORE -g gene_id -t transcript -p --countReadPairs"
    FC_EXONS_PARAMS = " -R CORE -g gene_id -J -p --countReadPairs"
    FC_TRANSCRIPTS_PARAMS = " -R CORE -g transcript_id -t exon -O -p --countReadPairs"
    FC_EXONBINS_PARAMS = " -R CORE -f -g exon_id -t exonic_part -O -p --countReadPairs"
    FC_EEJ_PARAMS = (
        " -R CORE -g junction_id -t eej -O -p --countReadPairs --fracOverlapFeature 0.9"
    )
    FC_EIJ_PARAMS = " -R CORE -g junction_id -t eij -O -p --countReadPairs"
else:
    FC_GENES_PARAMS = " -R CORE -g gene_id -t transcript"
    FC_EXONS_PARAMS = " -R CORE -g gene_id -J"
    FC_TRANSCRIPTS_PARAMS = " -R CORE -g transcript_id -t exon -O"
    FC_EXONBINS_PARAMS = " -R CORE -f -g exon_id -t exonic_part -O"
    FC_EEJ_PARAMS = " -R CORE -g junction_id -t eej -O --fracOverlapFeature 0.9"
    FC_EIJ_PARAMS = " -R CORE -g junction_id -t eij -O"


### Target rule input

if not (any(config["final_output"].values()) or config.get("mutpos", False)):
    raise ValueError("One of final_output values must be True!")


def get_other_output():
    target = []

    if config["final_output"].get("cB", False):
        target.append("results/cB/cB.csv.gz")

    if config["final_output"].get("cUP", False):
        target.append("results/cUP/cUP.csv.gz")

    if config["final_output"].get("arrow", False):
        target.append(
            expand(
                "results/arrow_dataset/sample={sample}/part-0.parquet",
                sample=SAMP_NAMES,
            )
        )

    if config.get("make_tracks", True):
        target.append(
            expand(
                "results/tracks/{sample}.{mut}.{id}.{strand}.tdf",
                sample=SAMP_NAMES,
                mut=Mutation_Types,
                id=[0, 1, 2, 3, 4, 5],
                strand=["pos", "min"],
            )
        )

    if config.get("features").get("threeputr", False):
        target.append(THREEPUTR_ANNOTATION)

    if not config["bam2bakr"]:
        target.append("results/multiqc/multiqc_report.html")

    if config["aligner"] == "star" and not config["bam2bakr"] and config["run_rsem"]:
        target.append(expand("results/rsem/{SID}.isoforms.results", SID=SAMP_NAMES))

    if config.get("mutpos", False):
        target.append("results/cB/mutpos.csv.gz")
        target.append("results/cB/mutpos_filtered.csv.gz")

    return target


STRAND = "R" if config["strandedness"] == "reverse" else "F"

MAKECB_THREADS = max(1, len(SAMP_NAMES) * 4)

CBINPUT = expand(
    "results/merge_features_and_muts/{sample}_cB.csv",
    sample=SAMP_NAMES,
)

CUPINPUT = expand(
    "results/merge_features_and_muts/{sample}_cUP.csv",
    sample=SAMP_NAMES,
)


##### MULTIQC INPUT ######

MULTIQC_INPUT = expand(
    "results/fastqc/{sample}_r{read}_fastqc.{ext}",
    sample=SAMP_NAMES,
    read=READS,
    ext=["html", "zip"],
)

if config.get("aligner") == "star":
    MULTIQC_INPUT.append(
        expand("results/align/{sample}-Log.final.out", sample=SAMP_NAMES)
    )
elif config.get("aligner") == "hisat3n":
    MULTIQC_INPUT.append(
        expand("results/align/{sample}_summary.txt", sample=SAMP_NAMES)
    )
