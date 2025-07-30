'''
    File name: mut_call.py
    Author: Martin Machyna
    Email: machyna@gmail.com
    Orcid: 0000-0002-3624-3472
    Wikidata: Q55137016
    Date created: 8/20/2021
    Date last modified: 8/30/2021
    Version: 1.0.0
    License: GPLv3
    Python Version: Python 2.7.18, 3.8.7+
    Packages: pysam 0.16.0.1
    Description: Script for detecting mutations from sequencing data .bam file
'''

import pysam
import csv
import argparse
import datetime

# Parse commandline arguments
parser = argparse.ArgumentParser(description='This is python implementation of TimeLapse mutation calling')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-b', '--bam', type=str, required=True, metavar = 'in_file.bam',
                    help='Bam file to process')
parser.add_argument('--mutType', default='TC', type=str,
                    help='Type of mutation to record (default: TC)')
parser.add_argument('--reads', default='PE', type=str, choices=['PE', 'SE'],
                    help='Type of mutation to record (default: PE)')
parser.add_argument('--minDist', default=-1, type=int, metavar = '<int>',
                    help='Base distance from read-end filter (default: -1)')
parser.add_argument('--minQual', default=40, type=int, metavar = '<int>',
                    help='Base minimal quality filter (default: 40)')
parser.add_argument('--tracks', action='store_true', # Automatically stored as default = FALSE
                    help='Generate files necessary for creating browser tracks')
parser.add_argument('--mutsRDS', action='store_true',
                    help='Generate _muts.rds like file with mutation frequency output')
parser.add_argument('--mutPos', action='store_true',
                    help='Generate _cU.rds like file and mutation position bedGraph files')
parser.add_argument('--SNPs', help='Path to snp.txt file')
parser.add_argument('--strandedness', default='F', type=str, choices=['F', 'R'],
                    help='Is first read forward or reverse orientation? F = forward is the default.')
args = parser.parse_args()

args.mutType = args.mutType.split(',')
args.base = [x[0] for x in args.mutType]        # base nucleotide: TC => T
inputName = args.bam.split('.bam')[0]           # name without .bam suffix

if args.strandedness == 'F':
    strand_check = True
else:
    strand_check = False



# Initialize variables
freq = {}               # dictionary for _muts.csv file with structure -> key='chrom:pos'  valuses=[trials, muts]
cU = {}
firstReadName = ''
muts = {'TA': 0, 'CA': 0, 'GA': 0, 'NA': 0, 'AT': 0, 'CT': 0, 'GT': 0, 'NT': 0, 'AC': 0, 'TC': 0, 'GC': 0, 'NC': 0, 'AG': 0, 'TG': 0, 'CG': 0, 'NG': 0, 'AN': 0, 'TN': 0, 'CN': 0, 'GN': 0, 'NN': 0}
DNAcode={'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n'}  # DNA code for comp and revcomp transformation
header = ['qname', 'nA', 'nC', 'nT', 'nG', 'rname', 'FR', 'sj', 'TA', 'CA', 'GA', 'NA', 'AT', 'CT', 'GT', 'NT', 'AC', 'TC', 'GC', 'NC', 'AG', 'TG', 'CG', 'NG', 'AN', 'TN', 'CN', 'GN', 'NN']

# For counting mutations at individual positions
if args.mutPos:
    header.extend(['gmutloc', 'tp'])


r_info = [''] + 4*[0] + 3*['']
dovetail = []
MDstore = {}


# Load SNPs for filtering
snp = {}
snpFile = open(args.SNPs, 'r')
for line in snpFile:
    line = line.strip().split(':')
    snp[line[2] + ':' + line[3]] = line[0] + ':' + line[1]

# Create files for tracks
if args.tracks:
    # Create output files names
    fileName = []
    for mt in args.mutType:
        for i in range(0,6):
            fileName.append('_'.join([inputName, mt, str(i), 'reads.txt']))

    # Open all files for writing
    fs = []
    for f in fileName:
        fs.append(open(f, 'w'))


# NO LONGER WRITING TO CU AS I ITERATE OVER BAM FILE
# TO SAVE ON DISK SPACE.
# # Initialize cU file
# if args.mutPos:
#     wr.writerow(['rname', 'gloc', 'GF', 'XF', 'ai', 'tp', 'trials', 'n'])
#
#     mycU = open(inputName + '_cU.csv', 'w', newline='')
#     cUwr = csv.writer(mycU)
#     cUwr.writerow(['rname', 'gloc', 'GF', 'XF', 'ai', 'tp', 'trials', 'n'])



# Set .bam file for reading
samfile = pysam.AlignmentFile(args.bam, 'rb')

print('Start: ' + str(datetime.datetime.now()))
for r in samfile:

    # Initialize + acquire info: First read only
    if firstReadName != r.query_name:
        muts={'TA': 0, 'CA': 0, 'GA': 0, 'NA': 0, 'AT': 0, 'CT': 0, 'GT': 0, 'NT': 0, 'AC': 0, 'TC': 0, 'GC': 0, 'NC': 0, 'AG': 0, 'TG': 0, 'CG': 0, 'NG': 0, 'AN': 0, 'TN': 0, 'CN': 0, 'GN': 0, 'NN': 0}
        r_info = [''] + 4*[0] + 3*['']
        dovetail = []
        MDstore = {}
        gmutloc = []
        tp = []

        r_info[0] = r.query_name            # Read name
        r_info[5] = r.reference_name        # Chromosome name


    # Gather alignmet information + Resolve dovetailing: Both reads
    if ('I' not in r.cigarstring) and ('D' not in r.cigarstring):       # Any read without insertions/deletions

        r_info[7] = str( r_info[7] == 'TRUE' or ('N' in r.cigarstring) ).upper()     # sj: splice junction

        if (r.is_paired and (r.is_read1 == (r.is_reverse == strand_check))) or (not r.is_paired and (r.is_reverse == strand_check)):        # If read is first_in_pair and on reverse strand -or- second_in_pair and on forward strand then make sequence complement
            r_info[6] = 'R'      # FR: forward or reverse read orientation
            MD = [[x[1], DNAcode[x[2]], min(x[0] - r.query_alignment_start + 1, r.query_alignment_length - (x[0] - r.query_alignment_start))] for x in r.get_aligned_pairs(matches_only = True, with_seq=True)]
            # Parse MD and Cigar strings, remove values that are softclipped
            # MD = [[gen_position, ref_base, base_readEnd_distance]]

            temp_qual = r.query_qualities
            r.query_sequence = ''.join([DNAcode[x] for x in r.query_sequence])
            r.query_qualities = temp_qual
        else:
            r_info[6] = 'F'
            MD = [[x[1], x[2], min(x[0] - r.query_alignment_start + 1, r.query_alignment_length - (x[0] - r.query_alignment_start))] for x in r.get_aligned_pairs(matches_only = True, with_seq=True)]



        if firstReadName != r.query_name:       # First read
            MDstore = {z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)}
            # store informatinon in dictionary of lists: {gen_position: [ref_base, read_base, qual, base_readEnd_distance]}
        else:                                   # Second read
            dovetail = list(set(MDstore.keys()) & set([x[0] for x in MD]))   # Identify genomic positions that are covered by both first and second read


            if len(dovetail) == 0:      # No dovetailing
                MDstore.update({z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)})

            else:                       # Dovetailing
                MD = {z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)}

                # MDstore.update({ pos:data for pos, data in MD.items() if pos in dovetail and MDstore[pos][2] < data[2] })   # Replace dovetail positions if better quality
                MDstore.update({ pos:data for pos, data in MD.items() if pos in dovetail and ((MDstore[pos][2] < data[2] and MDstore[pos][0].islower() and data[0].islower()) or (MDstore[pos][2] < data[2] and MDstore[pos][0].isupper() and data[0].isupper()) or (MDstore[pos][2] < data[2] and MDstore[pos][0].islower() and data[0].isupper() and MDstore[pos][2] + 33 < args.minQual) or (data[0].islower() and MDstore[pos][0].isupper() and data[2] + 33 > args.minQual)) })
                # This is a hack to simulate TimeLapse.R behaviour, but does not necessarily mean that it is a correct dovetail mutations handling
                # For dovetail bases: 1) If there is no mutation in 1st and in 2nd read => replace with higher quality 2nd read
                #                     2) If there is mutation in 1st and in 2nd read => replace with higher quality 2nd read
                #                     3) If there is mutation in 1nd but not in 2st read => replace only if 1st read quality is less than threshold and less than 2nd read
                #                     4) If there is mutation in 2nd read but not in 1st read => replace even with lower quality 2nd read as long as 2nd read quality is higher than threshold


                MDstore.update({ pos:data for pos, data in MD.items() if pos not in dovetail })                               # Append non dovetail sites



    # Collect data: Second read only or if in SE mode
    if (args.reads == 'SE' or firstReadName == r.query_name) and len(MDstore) > 0:
        refseq = [x[0].upper() for x in MDstore.values() if x[2] + 33 > args.minQual]  # Get reference sequence for readpair keeping only bases with given qaulity (Note: I think this should be also filtered for closeness to read end and presence of SNPs)
        # Count bases in reference sequence (soft clipped, dovetail-free)
        r_info[1] = refseq.count('A')       # nA
        r_info[2] = refseq.count('C')       # nC
        r_info[3] = refseq.count('T')       # nT
        r_info[4] = refseq.count('G')       # nG


        # Loop through every base of alignment and find mutations
        for pos, b in MDstore.items():

            # _muts.rds data
            if args.mutsRDS:
                if (b[0].upper() in args.base):
                    if (r.reference_name + ':' + str(pos)) not in freq:
                        freq[r.reference_name + ':' + str(pos)] = [1, 0]                   # Initialize counter for position
                    else:
                        freq[r.reference_name + ':' + str(pos)][0] += 1                    # Increment read coverage counter for given genomic position

                if b[0].islower() and ((b[0].upper() + b[1]) in args.mutType):
                    freq[r.reference_name + ':' + str(pos)][1] += 1                        # Increment mutation counter for given genomic position

            # cU.rds trial data
            if args.mutPos:
                if (b[0].upper() in args.base):
                    whichMut = [mut for mut in args.mutType if mut[0] == b[0].upper()]     # Find out which mutation types use this reference base e.g. T -> TC, TG, TA, TN
                    for mt in whichMut:
                        key = r.reference_name + ':' + str(pos) + ':' + mt
                        if key not in cU:
                            cU[key] = [1, 0]
                        else:
                            cU[key][0] += 1

            # _counts.rds data
            if b[0].islower() and (b[2] + 33 > args.minQual) and (b[3] > args.minDist) and (r.reference_name + ':' + str(pos + 1) not in snp):   # Find mutations marked as lowercase letters; apply quality filter; apply distance to read end filter; position is not a SNP
                muts[b[0].upper() + b[1]] += 1                                            # Increment the mutation counter for current readpair

                # mutPos bedGraph data + cU.rds n data
                if args.mutPos:
                    if (b[0].upper() + b[1]) in args.mutType:
                        key = r.reference_name + ':' + str(pos) + ':' + b[0].upper() + b[1]
                        cU[key][1] += 1

                        key = r.reference_name +  ':' + str(pos) + ':' + r_info[6] + ':' + b[0].upper() + b[1]

                        gmutloc.append(str(pos))            # Record position of muatation
                        tp.append(b[0].upper() + b[1])      # Record type of mutation




        # Write read info into _counts.csv
        r_info.extend( list(muts.values()) )
        if args.mutPos:
            r_info.extend( [ '|'.join(gmutloc), '|'.join(tp) ] )

        # Was doing this to save on RAM
        # Realized this inflated disk space usage tremendously
        # Will save here for posterity's sake
        # # Write to cU file
        # if args.mutPos:
        #
        #     for position, counts in cU.items():
        #         row = position.split(':')
        #         row[1] = int(row[1]) + 1                        # adjust position because we are 0-based
        #         row.extend(counts)
        #         cUwr.writerow(row)
        # 
        #     cU = {}

        # Save read name to track files
        if args.tracks:
            for index, mut in enumerate(args.mutType):
                for c in range(0, (5 if muts[mut] > 5 else muts[mut]) + 1):
                    fs[c + index * 6].write(r.query_name + '\n')

    # Save read name for next iteration
    firstReadName = r.query_name

print('end: ' + str(datetime.datetime.now()))


##### Close files ######
myfile.close()

### saving cU.rds file
if args.mutPos:
    with open(inputName + '_cU.csv', 'w', newline='') as myfile:
        wr = csv.writer(myfile)
        wr.writerow(['rname', 'gloc', 'tp', 'trials', 'n'])
        for position, counts in cU.items():
            row = position.split(':')
            row[1] = int(row[1]) + 1                        # adjust position because we are 0-based
            row.extend(counts)
            wr.writerow(row)

    del cU
print('cU: ' + str(datetime.datetime.now()))



if args.tracks:
    for f in fs:
        f.close()

