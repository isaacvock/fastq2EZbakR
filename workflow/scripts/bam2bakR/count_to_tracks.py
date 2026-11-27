'''
    File name: count_to_trcks.py
    Author: Martin Machyna
    Email: machyna@gmail.com
    Orcid: 0000-0002-3624-3472
    Wikidata: Q55137016
    Date created: 8/24/2021
    Date last modified: 8/24/2021
    Version: 1.0.0
    License: GPLv3
    Python Version: Python 2.7.18, 3.8.7+
    Description: Script for creating files contating names of reads that have given number of mutations. It supports any mutation types or their combination.
'''


import gzip
import csv
import argparse

from itertools import cycle
from itertools import chain

# Parse commandline arguments
parser = argparse.ArgumentParser(description='Script for creating files contating names of reads that bear given number of mutations')
parser.add_argument('-i', '--inputFile', type=str,
                    help='Name of gzip compressed csv file that will be processed')
parser.add_argument('-s', '--sample', type=str,
                    help='Sample name used for name of output file')
parser.add_argument('--mutType', default='TC', type=str,
                    help='Type of mutation to record (default: TC)')

args = parser.parse_args()
args.mutType = args.mutType.split(",")


# Create output files names
fileName = []
for b in args.mutType:
    for i in range(0,6):
        fileName.append("_".join([args.sample, b, str(i), "reads.txt"]))





with gzip.open(args.inputFile, mode="rt") as f:
    csv_reader = csv.reader(f)
    header = next(csv_reader)                                   # Read header
    mutIndex = [header.index(mut) for mut in args.mutType]      # Get column number for wanted mutation type(s)


    fs = []
    for f in fileName:
        fs.append(open(f, 'w'))                                 # Open all files for writing

    table = list(zip( range(0,len(fileName)), cycle(range(0,6)), chain(*[[x]*6 for x in mutIndex]) ))  # create list[[file_index, mutations_count, mutation_type_index], ...]


    try:
        for row in csv_reader:

            for z in table:
                if int( row[z[2]] ) >= z[1]: fs[z[0]].write(row[0] + '\n')   # if  mut_type_count >= count  then   fs[file_id].write(read_name + '\n')
    finally:
        for fh in fs:
            fh.close()
