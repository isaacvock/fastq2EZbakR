import glob
import os


### GENERAL HELPER FUNCTIONS/VARIABLES USED IN ALL CASES

# Sample names to help expanding lists of all bam files
# and to aid in defining wildcards
SAMP_NAMES = list(config['samples'].keys())

# Directory containing index; used in case of certain aligners
INDEX_DIR = config["indices"]

# Determine how many fastqs to look for
if config["PE"]:
    READS = [1, 2]
    READ_NAMES = ['r1', 'r2']
else:
    READS = [1]
    READ_NAMES = ['r1']

# Get input fastq files for first step
def get_input_fastqs(wildcards):
    fastq_path = config["samples"][wildcards.sample]
    fastq_files = sorted(glob.glob(f"{fastq_path}/*.fastq*"))
    return fastq_files

# Check if fastq files are gzipped
fastq_paths = config["samples"]

is_gz = False

for p in fastq_paths.values():

    fastqs = sorted(glob.glob(f"{p}/*.fastq*"))
    test_gz = any(path.endswith('.fastq.gz') for path in fastqs)
    is_gz = any([is_gz, test_gz])



### STAR HELPERS

# Trimmed fastq file paths, used as input for aligners

def get_fastq_r1(wildcards):
    if config["PE"]:

        return expand("results/trimmed/{SID}.1.fastq", SID = wildcards.sample)

    else:

        return expand("results/trimmed/{SID}.fastq", SID = wildcards.sample)

def get_fastq_r2(wildcards):
    if config["PE"]:

        return expand("results/trimmed/{SID}.2.fastq", SID = wildcards.sample)

    else:

        return ""


### Extra parameters passed to STAR alignment


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
sam_attributes = set(config["star_sam_tags"])

# Process --outSAMattributes
if "--outSAMattributes" in args:
    
    index = args.index("--outSAMattributes")

    # Assuming the attributes are space-separated and continuous after the flag
    next_flag_index = len(args)
    for i in range(index + 1, len(args)):
        if args[i].startswith('--'):
            next_flag_index = i
            break
    existing_attributes = set(args[index + 1:next_flag_index])
    combined_attributes = existing_attributes.union(sam_attributes)
    args[index + 1:next_flag_index] = list(combined_attributes)

else:
    args.extend(["--outSAMattributes"] + list(sam_attributes))

# Process --quantMode
quant_mode = set(["TranscriptomeSAM", "GeneCounts"])

if "--quantMode" in args:
    
    index = args.index("--quantMode")

    # Figure out what arguments were supplied to --quantMode
    next_flag_index = len(args)
    for i in range(index + 1, len(args)):
        if args[i].startswith('--'):
            next_flag_index = i
            break

    existing_attributes = set(args[index + 1:next_flag_index])

    # Add desired quantification modes
    combined_attributes = existing_attributes.union(quant_mode)
    args[index + 1:next_flag_index] = list(combined_attributes)

else:
    args.extend(["--quantMode"] + list(quant_mode))

# Force --outSAMtype to be BAM SortedByCoordinate
if "--outSAMType" in args:
    
    index = args.index("--outSAMType")

    # Replace the existing outSAMType values with the desired ones
    args[index + 1:index + 3] = ["BAM", "SortedByCoordinate"]

else:
    
    args.extend(["--outSAMType", "BAM", "SortedByCoordinate"])

# Force --sjdbGTFfile to be the provided annotation
if "--sjdbGTFfile" in args:
    
    index = args.index("--sjdbGTFfile")

    # Replace the existing outSAMType values with the desired ones
    args[index + 1:index + 3] = [str(config["annotation"])]

else:
    
    args.extend(["--sjdbGTFfile", str(config["annotation"])])



STAR_EXTRA = ' '.join(args)


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
CTL_NAMES = list(config['control_samples'])


# +s4U sample names
s4U_SAMPS = list(set(SAMP_NAMES) - set(CTL_NAMES))


# Number of -s4U control samples
nctl = len(CTL_NAMES)


# Number of +s4U samples
nsamps = len(SAMP_NAMES)


# Name for genome index
def get_index_name():
    genome = config["genome_fasta"]
    index = str(genome) + ".fai"
    return index


# PE or SE? 

if config["PE"]:

    FORMAT = "PE"

else:

    FORMAT = "SE"


# Calculate scale factors for tracks?
NORMALIZE = config['normalize']



# bam2EZbakR bam file fetching
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]


# pnew estimates for RSEM+
def get_pnew(wildcards):
    return config["pnews"][wildcards.sample]

# pold estimates for RSEM+
def get_pold(wildcards):
    return config["polds"][wildcards.sample]

