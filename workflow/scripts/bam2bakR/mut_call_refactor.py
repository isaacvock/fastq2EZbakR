#!/usr/bin/env python3
import pysam
import csv
import argparse
import datetime
import gzip                         # (3) write gzip directly
from collections import Counter     # (4) fast nucleotide tallies

# ------------------------------------------------------------------
# Parse command‑line arguments – 1 new option (--threads)
# ------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description='Python implementation of TimeLapse mutation calling')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-b', '--bam', type=str, required=True,
                           metavar='in_file.bam', help='Bam file to process')
parser.add_argument("-o", "--output", required=True,
                   help="gzip‑compressed per‑read CSV")
parser.add_argument('--threads', default=1, type=int, metavar='<int>',    # NEW
                    help='BGZF decompression threads (default: 1)')
parser.add_argument('--reads', default='PE', type=str, choices=['PE', 'SE'],
                    help='Type of mutation to record (default: PE)')
parser.add_argument('--minDist', default=-1, type=int, metavar='<int>',
                    help='Base distance from read‑end filter (default: -1)')
parser.add_argument('--minQual', default=40, type=int, metavar='<int>',
                    help='Base minimal quality filter (default: 40)')
parser.add_argument('--SNPs', help='Path to snp.txt file')
parser.add_argument('--strandedness', default='F', type=str,
                    choices=['F', 'R'],
                    help='Is first read forward (F) or reverse (R) orientation?')
parser.add_argument('--mutPos', action='store_true',                    # kept
                    help='Output per‑mutation position columns')
args = parser.parse_args()

strand_check = True if args.strandedness == 'F' else False

# ------------------------------------------------------------------
# (2) translation table for fast complement
# ------------------------------------------------------------------
DNA_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")   # replaces dict look‑ups
DNAcode    = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N',
              'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n'}  # kept for MD

# ------------------------------------------------------------------
# Initialize variables – identical to original
# ------------------------------------------------------------------
freq = {}
cU = {}
firstReadName = ''
muts = {'TA': 0, 'CA': 0, 'GA': 0, 'NA': 0, 'AT': 0, 'CT': 0, 'GT': 0, 'NT': 0,
        'AC': 0, 'TC': 0, 'GC': 0, 'NC': 0, 'AG': 0, 'TG': 0, 'CG': 0, 'NG': 0,
        'AN': 0, 'TN': 0, 'CN': 0, 'GN': 0, 'NN': 0}
header = ['qname', 'nA', 'nC', 'nT', 'nG', 'rname', 'FR', 'sj'] + list(muts.keys())

if args.mutPos:
    header.extend(['gmutloc', 'tp'])

r_info   = [''] + 4*[0] + 3*['']
dovetail = []
MDstore  = {}

# ------------------------------------------------------------------
# Load SNP table
# ------------------------------------------------------------------
snp = {}
snpFile = open(args.SNPs, 'r')
for line in snpFile:
    line = line.strip().split(':')
    snp[line[2] + ':' + line[3]] = line[0] + ':' + line[1]

# ------------------------------------------------------------------
# (3) Write gzip directly – output file now ends in .csv.gz
# ------------------------------------------------------------------
with gzip.open(args.output, "wt", newline="") as myfile:
    wr = csv.writer(myfile)
    wr.writerow(header)

    # ------------------------------------------------------------------
    # (1) open BAM with multithreaded BGZF decompression
    #     (6) allows "-" for stdin if you want to pipe samtools view
    # ------------------------------------------------------------------
    bam_source = "-" if args.bam == "-" else args.bam
    samfile = pysam.AlignmentFile(bam_source, 'rb', threads=args.threads)

    print('Start:', datetime.datetime.now())
    for r in samfile:

        # -------- Initialize per‑read‑pair vars (unchanged) --------
        if firstReadName != r.query_name:
            muts = dict.fromkeys(muts, 0)   # reset counts
            r_info = [''] + 4*[0] + 3*['']
            dovetail = []
            MDstore = {}
            gmutloc = []
            tp = []
            r_info[0] = r.query_name
            r_info[5] = r.reference_name

        # -------- Gather alignment info + dovetail handling --------
        if ('I' not in r.cigarstring) and ('D' not in r.cigarstring):

            r_info[7] = str(r_info[7] == 'TRUE' or ('N' in r.cigarstring)).upper()

            do_complement = ((r.is_paired and (r.is_read1 == (r.is_reverse == strand_check))) or
                             (not r.is_paired and (r.is_reverse == strand_check)))

            if do_complement:
                r_info[6] = 'R'
                temp_qual = r.query_qualities
                r.query_sequence = r.query_sequence.translate(DNA_TRANS)      # (2)
                r.query_qualities = temp_qual
                MD = [[x[1], DNAcode[x[2]], min(x[0]-r.query_alignment_start+1,
                                                r.query_alignment_length-(x[0]-r.query_alignment_start))]
                      for x in r.get_aligned_pairs(matches_only=True, with_seq=True)]
            else:
                r_info[6] = 'F'
                MD = [[x[1], x[2], min(x[0]-r.query_alignment_start+1,
                                        r.query_alignment_length-(x[0]-r.query_alignment_start))]
                      for x in r.get_aligned_pairs(matches_only=True, with_seq=True)]

            # ---- dovetail logic remains byte‑for‑byte identical ----
            if firstReadName != r.query_name:
                MDstore = {z[0][0]: [z[0][1], z[1], z[2], z[0][2]]
                           for z in zip(MD, r.query_alignment_sequence,
                                        r.query_alignment_qualities)}
            else:
                dovetail = list(set(MDstore.keys()) &
                                set([x[0] for x in MD]))
                if len(dovetail) == 0:
                    MDstore.update({z[0][0]: [z[0][1], z[1], z[2], z[0][2]]
                                    for z in zip(MD, r.query_alignment_sequence,
                                                 r.query_alignment_qualities)})
                else:
                    MD = {z[0][0]: [z[0][1], z[1], z[2], z[0][2]]
                          for z in zip(MD, r.query_alignment_sequence,
                                       r.query_alignment_qualities)}
                    MDstore.update({pos: data for pos, data in MD.items()
                                    if pos in dovetail and
                                    ((MDstore[pos][2] < data[2] and
                                      MDstore[pos][0].islower() and data[0].islower()) or
                                     (MDstore[pos][2] < data[2] and
                                      MDstore[pos][0].isupper() and data[0].isupper()) or
                                     (MDstore[pos][2] < data[2] and
                                      MDstore[pos][0].islower() and data[0].isupper() and
                                      MDstore[pos][2] + 33 < args.minQual) or
                                     (data[0].islower() and MDstore[pos][0].isupper() and
                                      data[2] + 33 > args.minQual))})
                    MDstore.update({pos: data for pos, data in MD.items()
                                    if pos not in dovetail})

        # -------- Collect data: second read only (unchanged) -------
        if (args.reads == 'SE' or firstReadName == r.query_name) and MDstore:
            refseq = [x[0].upper() for x in MDstore.values()
                      if x[2] + 33 > args.minQual]

            # (4) Counter replaces 4 × list.count
            cnt = Counter(refseq)
            r_info[1] = cnt['A']
            r_info[2] = cnt['C']
            r_info[3] = cnt['T']
            r_info[4] = cnt['G']

            for pos, b in MDstore.items():
                if (b[0].islower() and (b[2] + 33 > args.minQual) and
                    (b[3] > args.minDist) and
                    (f"{r.reference_name}:{pos+1}" not in snp)):
                    muts[b[0].upper() + b[1]] += 1

            r_info.extend(list(muts.values()))
            if args.mutPos:
                r_info.extend(['|'.join(gmutloc), '|'.join(tp)])
            wr.writerow(r_info)

        firstReadName = r.query_name

    print('End:', datetime.datetime.now())
    samfile.close()   # automatically closed on context‑exit too
