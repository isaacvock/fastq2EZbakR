import glob
import os

### Figure out if bam files or fastq files are provided
config["bam2bakr"] = all(value.endswith(".bam") for value in config["samples"].values())


### If assigning reads to junctions, need to keep jI and jM tags for now
if config["features"]["junctions"]:
    config["remove_tags"] = False


### GENERAL HELPER FUNCTIONS/VARIABLES USED IN ALL CASES

# Sample names to help expanding lists of all bam files
# and to aid in defining wildcards
if config["download_fastqs"]:
    SAMP_NAMES = config["sra_accessions"]
else:
    SAMP_NAMES = list(config["samples"].keys())

NUM_SAMPS = len(SAMP_NAMES)


# Directory containing index; used in case of certain aligners
INDEX_DIR = config["indices"]


# Make life easier for users and catch if they add a '/' at the end of their path
# to alignment indices. If so, remove it to avoid double '/'
if config["indices"].endswith("/"):
    INDEX_PATH = str(config["indices"])
    INDEX_PATH = INDEX_PATH[:-1]
else:
    INDEX_PATH = str(config["indices"])


# Determine how many fastqs to look for
if config["PE"]:
    READS = [1, 2]
    READ_NAMES = ["r1", "r2"]
else:
    READS = [1]
    READ_NAMES = ["r1"]


# Define the function to get input FASTQ files
# Suggested in this issue by ld32: https://github.com/isaacvock/fastq2EZbakR/issues/23
def get_input_fastqs(wildcards):
    print(f"Resolving input for sample: {wildcards.sample}")

    if config.get("download_fastqs", True):
        print("Download fastqs: Yes")
        if config.get("PE", False):
            SRA_READS = ["_1.fastq", "_2.fastq"]
            print(f"Sample {wildcards.sample} is paired-end.")
        else:
            SRA_READS = [".fastq"]
            print(f"Sample {wildcards.sample} is single-end.")

        paths = [
            f"results/download_fastq/{wildcards.sample} {read}" for read in SRA_READS
        ]
        print(f"Paths: {sorted(paths)}")
        return sorted(paths)

    else:
        print("Download fastqs: No")
        fastq_path = config["samples"].get(wildcards.sample, None)
        if fastq_path is None:
            raise ValueError(f"No path found for sample {wildcards.sample}")
        print(f"path is: {fastq_path}")
        current_directory = os.getcwd()
        print("Current Working Directory:", current_directory)
        all_files = os.listdir(fastq_path)
        print(f"All files in directory: {all_files}")

        fastq_files = []

        # List all entries in the directory and store full paths in the list
        for entry in os.listdir(fastq_path):
            full_path = os.path.join(fastq_path, entry)
            fastq_files.append(full_path)

        # Filter files to only include those that end with .fastq or .fastq.gz
        # fastq_files = [f"data/fastq/{wildcards.sample}/{wildcards.sample}.{n}.fastq.gz" for n in ("1", "2")] #sorted([os.path.join(directory, f) for f in all_files if f.endswith('.fastq') or f.endswith('.fastq.gz')])
        # fastq_files = sorted(glob.glob(f"{fastq_path}/*.fastq*"))
        print(f"Found files: {sorted(fastq_files)}")
        return sorted(fastq_files)


# Check if fastq files are gzipped
fastq_paths = config["samples"]


if not config["download_fastqs"]:
    is_gz = False

    for p in fastq_paths.values():
        fastqs = sorted(set(glob.glob(f"{p}/*.fastq*") + glob.glob(f"{p}/*.fq*")))
        test_gz = any(path.endswith(".fastq.gz") for path in fastqs)
        is_gz = any([is_gz, test_gz])


# Columns in final cB
keepcols = ["sample", "sj", "rname"]

if config["features"]["genes"]:
    keepcols.append("GF")

if config["features"]["exons"]:
    keepcols.append("XF")

if config["features"]["transcripts"]:
    keepcols.append("transcripts")

if config["features"]["exonic_bins"]:
    keepcols.append("exon_bin")

if config["features"]["tec"]:
    keepcols.append("TEC")

if config["features"]["junctions"]:
    keepcols.append("junction_start")
    keepcols.append("junction_end")

if config.get("features").get("threeputr", False):
    keepcols.append("threepUTR")

if config["features"]["eij"]:
    keepcols.append("ei_junction_id")

if config["features"]["eej"]:
    keepcols.append("ee_junction_id")


# Get mutation types to track
MutTypes = config["mut_tracks"]
Mutation_Types = MutTypes.split(",")
Nucleotide_Types = ["n{}".format(muttype[0]) for muttype in Mutation_Types]

# Columns I want to summarise by
init_keepcols = keepcols.copy()
cols_to_search = keepcols.copy()
cols_to_search.extend(Mutation_Types)
cols_to_search.extend(Nucleotide_Types)

keepcols = ",".join(keepcols)


### FASTQC HELPERS


def get_fastqc_read(wildcards):
    if config["skip_trimming"]:
        if is_gz:
            return expand(
                "results/unzipped/{SID}.{READ}.fastq",
                SID=wildcards.sample,
                READ=wildcards.read,
            )

        else:
            fastq_path = config["samples"][wildcards.sample]
            fastq_files = sorted(
                set(
                    glob.glob(f"{fastq_path}/*.fastq*")
                    + glob.glob(f"{fastq_path}/*.fq*")
                )
            )
            readID = int(wildcards.read) - 1
            return fastq_files[readID]

    else:
        return expand(
            "results/trimmed/{SID}.{READ}.fastq",
            SID=wildcards.sample,
            READ=wildcards.read,
        )


### STAR HELPERS

# Trimmed fastq file paths, used as input for aligners


def get_fastq_r1(wildcards):
    if config["skip_trimming"]:
        if is_gz:
            return expand("results/unzipped/{SID}.1.fastq", SID=wildcards.sample)

        else:
            fastq_path = config["samples"][wildcards.sample]
            fastq_files = sorted(
                set(
                    glob.glob(f"{fastq_path}/*.fastq*")
                    + glob.glob(f"{fastq_path}/*.fq*")
                )
            )
            return fastq_files[0]

    elif config["PE"]:
        return expand("results/trimmed/{SID}.1.fastq", SID=wildcards.sample)

    else:
        return expand("results/trimmed/{SID}.1.fastq", SID=wildcards.sample)


def get_fastq_r2(wildcards):
    if config["skip_trimming"]:
        if is_gz:
            return expand("results/unzipped/{SID}.2.fastq", SID=wildcards.sample)

        else:
            fastq_path = config["samples"][wildcards.sample]
            fastq_files = sorted(
                set(
                    glob.glob(f"{fastq_path}/*.fastq*")
                    + glob.glob(f"{fastq_path}/*.fq*")
                )
            )
            return fastq_files[1]

    elif config["PE"]:
        return expand("results/trimmed/{SID}.2.fastq", SID=wildcards.sample)

    else:
        return expand("results/trimmed/{SID}.2.fastq", SID=wildcards.sample)


# Determine whether or not indexing needs to be done
star_indices = glob.glob(f"{INDEX_PATH}/genomeParameters.txt")


if config["aligner"] == "star":
    if len(star_indices) == 0:
        make_index = True
    else:
        make_index = False


### HISAT2 HELPERS


def get_hisat2_reads(wildcards, READS=READS):
    if config["skip_trimming"]:
        if is_gz:
            return expand(
                "results/unzipped/{SID}.{read}.fastq", SID=wildcards.sample, read=READS
            )

        else:
            fastq_path = config["samples"][wildcards.sample]
            fastq_files = sorted(
                set(
                    glob.glob(f"{fastq_path}/*.fastq*")
                    + glob.glob(f"{fastq_path}/*.fq*")
                )
            )
            return fastq_files

    else:
        return expand(
            "results/trimmed/{SID}.{read}.fastq", SID=wildcards.sample, read=READS
        )


# Determine whether or not indexing needs to be done
hisat2_indices = glob.glob(f"{INDEX_PATH}/*.ht2*")

if config["aligner"] == "hisat2":
    if len(hisat2_indices) == 0:
        make_index = True
    else:
        make_index = False


### Extra parameters passed to STAR alignment


# Annotation to use for alignment and quantification
if config["quantify_premRNA"]:
    AandQ_ANNOTATION = "results/modify_annotation/modified_annotation.gtf"

else:
    AandQ_ANNOTATION = config["annotation"]

STAR_PARAMS = str(config["star_align_params"])

## Add necessary tags and change other settings as necessary

args = STAR_PARAMS.split()

# Tags to always include
# Technically, NH, HI, AS, and nM are added by default
# NH = number of places the read aligns
# HI = Numerical index that is useful to separate multi-mapping reads
# NM = Edit distance between read and reference
# AS = Alignment score
# MD = Used for mutation counting
# nM = Number of mismatches

if config["features"]["junctions"]:
    tags = list(config["star_sam_tags"])

    if "jI" not in tags:
        tags = tags + ["jI"]

    if "jM" not in tags:
        tags = tags + ["jM"]

    sam_attributes = set(tags)

else:
    sam_attributes = set(config["star_sam_tags"])

# Process --outSAMattributes
if "--outSAMattributes" in args:
    index = args.index("--outSAMattributes")

    # Assuming the attributes are space-separated and continuous after the flag
    next_flag_index = len(args)
    for i in range(index + 1, len(args)):
        if args[i].startswith("--"):
            next_flag_index = i
            break
    existing_attributes = set(args[index + 1 : next_flag_index])
    combined_attributes = existing_attributes.union(sam_attributes)
    args[index + 1 : next_flag_index] = list(combined_attributes)

else:
    args.extend(["--outSAMattributes"] + list(sam_attributes))

# Process --quantMode
quant_mode = set(["TranscriptomeSAM", "GeneCounts"])

if "--quantMode" in args:
    index = args.index("--quantMode")

    # Figure out what arguments were supplied to --quantMode
    next_flag_index = len(args)
    for i in range(index + 1, len(args)):
        if args[i].startswith("--"):
            next_flag_index = i
            break

    existing_attributes = set(args[index + 1 : next_flag_index])

    # Add desired quantification modes
    combined_attributes = existing_attributes.union(quant_mode)
    args[index + 1 : next_flag_index] = list(combined_attributes)

else:
    args.extend(["--quantMode"] + list(quant_mode))

# Force --outSAMtype to be BAM SortedByCoordinate
if "--outSAMtype" in args:
    index = args.index("--outSAMtype")

    # Replace the existing outSAMtype values with the desired ones
    args[index + 1 : index + 3] = ["BAM", "SortedByCoordinate"]

else:
    args.extend(["--outSAMtype", "BAM", "SortedByCoordinate"])

# Force --sjdbGTFfile to be the provided annotation
if "--sjdbGTFfile" in args:
    index = args.index("--sjdbGTFfile")

    # Replace the existing outSAMtype values with the desired ones
    args[index + 1 : index + 3] = [str(config["annotation"])]

else:
    args.extend(["--sjdbGTFfile", str(config["annotation"])])


STAR_EXTRA = " ".join(args)


### HISAT2 HELPERS

# Figure out what to pass to --rna-strandedness
if config["strandedness"] == "no":
    HISAT2_STRANDEDNESS = ""

elif config["strandedness"] == "reverse":
    if config["PE"]:
        HISAT2_STRANDEDNESS = "--rna-strandness RF"

    else:
        HISAT2_STRANDEDNESS = "--rna-strandness R"


else:
    if config["PE"]:
        HISAT2_STRANDEDNESS = "--rna-strandness FR"

    else:
        HISAT2_STRANDEDNESS = "--rna-strandness F"

# Index base name
if config["hisat2_index_base"]:
    HISAT2_BASE = "{}/{}".format(INDEX_PATH, config["hisat2_index_base"])

else:
    HISAT2_BASE = "{}/{}".format(INDEX_PATH, "hisat2_index")


### BAM2EZBAKR PARAMETERS


# -s4U control sample names
CTL_NAMES = list(config["control_samples"])


# +s4U sample names
s4U_SAMPS = list(set(SAMP_NAMES) - set(CTL_NAMES))


# Number of -s4U control samples
nctl = len(CTL_NAMES)


# Number of +s4U samples
nsamps = len(SAMP_NAMES)


# Name for genome index
def get_index_name():
    genome = config["genome"]
    index = str(genome) + ".fai"
    return index


# PE or SE?

if config["PE"]:
    FORMAT = "PE"

else:
    FORMAT = "SE"


# Calculate scale factors for tracks?
NORMALIZE = config["normalize"]


# bam2EZbakR bam file fetching
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


# pnew estimates for RSEM+
def get_pnew(wildcards):
    return config["pnews"][wildcards.sample]


# pold estimates for RSEM+
def get_pold(wildcards):
    return config["polds"][wildcards.sample]


### TRACK NORMALIZATION HELPERS

# What is input for normalization?
if config.get("use_exons_only", True):
    NORMALIZATION_INPUT = expand(
        "results/featurecounts_exons/{sample}.featureCounts", sample=SAMP_NAMES
    )

else:
    NORMALIZATION_INPUT = expand(
        "results/featurecounts_genes/{sample}.featureCounts", sample=SAMP_NAMES
    )

# Are exons being used?
if config.get("use_exons_only", True):
    NORMALIZATION_EXTRA = ""

else:
    NORMALIZATION_EXTRA = "--use_genes"


### FEATURECOUNTS HELPERS


## All of the files to merge
def get_merge_input(wildcards):
    MERGE_INPUT = []

    MERGE_INPUT.extend(
        expand("results/counts/{SID}_counts.csv.gz", SID=wildcards.sample)
    )

    if config.get("features").get("genes", True):
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_genes/{SID}.s.bam.featureCounts",
                SID=wildcards.sample,
            )
        )

    if config.get("features").get("exons", True):
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_exons/{SID}.s.bam.featureCounts",
                SID=wildcards.sample,
            )
        )

    if config.get("features").get("transcripts", False):
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_transcripts/{SID}.s.bam.featureCounts",
                SID=wildcards.sample,
            )
        )

    if config.get("features").get("exonic_bins", False):
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_exonbins/{SID}.s.bam.featureCounts",
                SID=wildcards.sample,
            )
        )

    if config.get("features").get("tec", False):
        MERGE_INPUT.extend(
            expand("results/read_to_transcripts/{SID}.csv", SID=wildcards.sample)
        )

    if config.get("features").get("junctions", False):
        MERGE_INPUT.extend(
            expand("results/read_to_junctions/{SID}.csv.gz", SID=wildcards.sample)
        )

    if config.get("features").get("threeputr", False):
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_3utr/{SID}.s.bam.featureCounts",
                SID=wildcards.sample,
            )
        )

    if config.get("features").get("eej", False):
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_eej/{SID}.s.bam.featureCounts",
                SID=wildcards.sample,
            )
        )

    if config.get("features").get("eij", False):
        MERGE_INPUT.extend(
            expand(
                "results/featurecounts_eij/{SID}.s.bam.featureCounts",
                SID=wildcards.sample,
            )
        )

    return MERGE_INPUT


# Get strandedness parameter
if config["strandedness"] == "reverse":
    FC_STRAND = 2

elif config["strandedness"] == "yes" or config["strandedness"] == "forward":
    FC_STRAND = 1

else:
    FC_STRAND = 0


## Get extra parameters for gene calling

if config["PE"]:
    FC_GENES_PARAMS = " -R CORE -g gene_id -t transcript -p --countReadPairs"

else:
    FC_GENES_PARAMS = " -R CORE -g gene_id -t transcript"


## Get extra parameters for exon calling

if config["PE"]:
    FC_EXONS_PARAMS = " -R CORE -g gene_id -J -p --countReadPairs"

else:
    FC_EXONS_PARAMS = " -R CORE -g gene_id -J"


## Get extra parameters for transcript calling

if config["PE"]:
    FC_TRANSCRIPTS_PARAMS = " -R CORE -g transcript_id -t exon -O -p --countReadPairs"

else:
    FC_TRANSCRIPTS_PARAMS = " -R CORE -g transcript_id -t exon -O"


## Get extra parameters for exon bin calling


if config["PE"]:
    FC_EXONBINS_PARAMS = " -R CORE -f -g exon_id -t exonic_part -O -p --countReadPairs"

else:
    FC_EXONBINS_PARAMS = " -R CORE -f -g exon_id -t exonic_part -O"


## Get extra parameters for exon-exon junction calling

if config["PE"]:
    FC_EEJ_PARAMS = (
        " -R CORE -g junction_id -t eej -O -p --countReadPairs --fracOverlapFeature 0.9"
    )

else:
    FC_EEJ_PARAMS = " -R CORE -g junction_id -t eej -O --fracOverlapFeature 0.9"


## Get extra parameters for exon-intron junction calling

if config["PE"]:
    FC_EIJ_PARAMS = " -R CORE -g junction_id -t eij -O -p --countReadPairs"

else:
    FC_EIJ_PARAMS = " -R CORE -g junction_id -t eij -O"


### Target rule input

# Need to specify one of final_output to be True
if not (any(config["final_output"].values()) or config.get("mutpos", False)):
    raise ValueError("One of final_output values must be True!")

# lowRAM can only work with a single type of output
if config["lowRAM"]:
    if config["final_output"]["cUP"]:
        raise ValueError("lowRAM = True is not currently compatible with cUP output!")


def get_other_output():
    target = []

    if config["final_output"]["cB"]:
        target.append("results/cB/cB.csv.gz")

    if config["final_output"]["cUP"]:
        target.append("results/cUP/cUP.csv.gz")

    if config["final_output"]["arrow"]:
        if config["lowRAM"]:
            target.append(
                expand(
                    "results/arrow_dataset/sample={sample}/part-0.csv",
                    sample=SAMP_NAMES,
                )
            )

        else:
            target.append(
                expand(
                    "results/arrow_dataset/sample={sample}/part-0.parquet",
                    sample=SAMP_NAMES,
                )
            )

    # Tracks
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

    # fastQC output
    if not config["bam2bakr"]:
        target.append("results/multiqc/multiqc_report.html")

    if config["aligner"] == "star" and not config["bam2bakr"] and config["run_rsem"]:
        target.append(expand("results/rsem/{SID}.isoforms.results", SID=SAMP_NAMES))

    if config["mutpos"]:
        target.append("results/cB/mutpos.csv.gz")

        target.append("results/cB/mutpos_filtered.csv.gz")

    # if config["use_salmon"]:

    #     target.append(expand("reul"))

    return target


### Strandedness parameter to pass to mutation counting
if config["strandedness"] == "reverse":
    STRAND = "R"

else:
    STRAND = "F"


### Optimizing cB creation multithreading

MAKECB_THREADS = len(SAMP_NAMES) * 4


### If keeping RAM usage low, need to determine which columns of merged files to support

# Columns as they appear in mutation counting csv output
colnames = [
    "qname",
    "nA",
    "nC",
    "nT",
    "nG",
    "rname",
    "FR",
    "sj",
    "TA",
    "CA",
    "GA",
    "NA",
    "AT",
    "CT",
    "GT",
    "NT",
    "AC",
    "TC",
    "GC",
    "NC",
    "AG",
    "TG",
    "CG",
    "NG",
    "AN",
    "TN",
    "CN",
    "GN",
    "NN",
]

if config["features"]["genes"]:
    colnames.append("GF")

if config["features"]["exons"]:
    colnames.append("XF")

if config["features"]["transcripts"]:
    colnames.append("transcripts")


if config["features"]["exonic_bins"]:
    colnames.append("exon_bin")

if config["features"]["tec"]:
    colnames.append("TEC")


if config["features"]["junctions"]:
    colnames.append("junction_start")
    colnames.append("junction_end")

if config.get("features").get("threeputr", False):
    colnames.append("threepUTR")


if config["features"]["eej"]:
    colnames.append("ee_junction_id")


if config["features"]["eij"]:
    colnames.append("ei_junction_id")


# Get indices of columns that I need to sort in order to summarise by
cols_to_sort = [index for index, item in enumerate(colnames) if item in cols_to_search]
cols_to_sort = sorted(cols_to_sort)

numericsort_cols = Mutation_Types.copy()
numericsort_cols.extend(Nucleotide_Types)
numeric_sort_columns = [
    index for index, item in enumerate(colnames) if item in numericsort_cols
]

# Need to add 1 because linux sort is 1-indexed
cols_to_sort = [x + 1 for x in cols_to_sort]
numeric_sort_columns = [x + 1 for x in numeric_sort_columns]

key_args = " ".join(
    [
        "-k{0},{0}{1}".format(col, "n" if col in numeric_sort_columns else "V")
        for col in cols_to_sort
    ]
)

SORTPARAMS = key_args

# Columns that will be summarized over, properly ordered
cols_to_summarize = [col for col in colnames if col in cols_to_search]
COLS_TO_SUM = ",".join(cols_to_summarize)


### Input for makecB

if config["lowRAM"]:
    if config["final_output"]["arrow"]:
        CBINPUT = expand(
            "results/arrow_dataset/sample={sample}/part-0.csv",
            sample=SAMP_NAMES,
        )

    else:
        CBINPUT = expand(
            "results/lowram_summarise/{sample}.csv",
            sample=SAMP_NAMES,
        )


else:
    CBINPUT = expand(
        "results/merge_features_and_muts/{sample}_cB.csv",
        sample=SAMP_NAMES,
    )


### Input for makecUP

if config["lowRAM"]:
    CUPINPUT = expand(
        "results/lowram_summarise/{sample}_cB.csv",
        sample=SAMP_NAMES,
    )


else:
    CUPINPUT = expand(
        "results/merge_features_and_muts/{sample}_cUP.csv",
        sample=SAMP_NAMES,
    )

### RSEM plus input


def get_rsemp_input(wildcards):
    RSEMP_INPUT = []

    if config["lowRAM"]:
        RSEMP_INPUT.extend(
            expand(
                "results/lowram_merge_features_and_counts/{SID}.csv",
                SID=wildcards.sample,
            )
        )

    else:
        RSEMP_INPUT.extend(
            expand(
                "results/merge_features_and_muts/{SID}_counts.csv.gz",
                SID=wildcards.sample,
            )
        )

    return RSEMP_INPUT


### Files to merge if lowRAM


## All of the files to merge
def get_lowram_merge_input(wildcards):
    MERGE_INPUT = []

    MERGE_INPUT.extend(
        expand("results/sort_mutcounts_by_qname/{SID}_counts.csv", SID=wildcards.sample)
    )

    if config["features"]["genes"]:
        MERGE_INPUT.extend(
            expand(
                "results/sort_fcgene_by_qname/{SID}.featureCounts", SID=wildcards.sample
            )
        )

    if config["features"]["exons"]:
        MERGE_INPUT.extend(
            expand(
                "results/sort_fcexon_by_qname/{SID}.featureCounts", SID=wildcards.sample
            )
        )

    if config["features"]["transcripts"]:
        MERGE_INPUT.extend(
            expand(
                "results/sort_fctranscript_by_qname/{SID}.featureCounts",
                SID=wildcards.sample,
            )
        )

    if config["features"]["exonic_bins"]:
        MERGE_INPUT.extend(
            expand(
                "results/sort_fcexonbin_by_qname/{SID}.featureCounts",
                SID=wildcards.sample,
            )
        )

    if config["features"]["tec"]:
        MERGE_INPUT.extend(
            expand(
                "results/sort_bamtranscript_by_qname/{SID}.csv", SID=wildcards.sample
            )
        )

    if config["features"]["junctions"]:
        MERGE_INPUT.extend(
            expand(
                "results/sort_junction_by_qname/{SID}_counts.csv", SID=wildcards.sample
            )
        )

    if config["features"]["eej"]:
        MERGE_INPUT.extend(
            expand(
                "results/sort_fcee_by_qname/{SID}.featureCounts", SID=wildcards.sample
            )
        )

    if config["features"]["eij"]:
        MERGE_INPUT.extend(
            expand(
                "results/sort_fcei_by_qname/{SID}.featureCounts", SID=wildcards.sample
            )
        )

    return MERGE_INPUT


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


##### 3'-UTR FEATURECOUNTS PARAMETER #####

if config["PE"]:
    if config.get("strandedness") == "reverse":
        FC_3UTR_PARAMS = "-R CORE -t 3UTR -O -p --read2pos 5"
    else:
        FC_3UTR_PARAMS = "-R CORE -t 3UTR -O -p --read2pos 3"
else:
    if config.get("strandedness") == "reverse":
        FC_3UTR_PARAMS = "-R CORE -t 3UTR -O --read2pos 5"
    else:
        FC_3UTR_PARAMS = "-R CORE -t 3UTR -O --read2pos 3"


if config.get("call_threeputrs", False):
    THREEPUTR_ANNOTATION = config.get("annotation")

else:
    THREEPUTR_ANNOTATION = "annotations/threepUTR_annotation.gtf"
